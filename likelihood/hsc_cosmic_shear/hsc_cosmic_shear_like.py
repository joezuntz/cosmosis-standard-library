import sys
import os
# add directory above this one to path
this_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(this_dir)
sacc_dir = os.path.join(parent_dir, "sacc")
sys.path.append(sacc_dir)
import sacc_like
from cosmosis.datablock import BlockError
import numpy as np


class HSCLikelihood(sacc_like.SaccClLikelihood):
    # This is a sub-class of the main SACC likelihood that adds
    # subtraction of a PSF-related template from the theory vector
    like_name = "hsc"

    def __init__(self, options):
        super().__init__(options)
        self.psf_ell_bins, self.psf_template, self.psf_transform_matrix, self.psf_means = self.load_psf_template(options)

    def load_psf_template(self, options):
        filename = options.get_string("psf_file", "")
        psf_ell_bins = np.load(filename)['arr_0']
        psf_cl_arr = np.load(filename)['arr_1']
        filename2 = options.get_string("psf_transformation_file", "")
        psf_transform_matrix = np.load(filename2)['arr_0']
        psf_means = np.load(filename2)['arr_1']
        return psf_ell_bins, psf_cl_arr, psf_transform_matrix, psf_means

    def compute_psf_template(self, block):
        alpha2 = block["psf_parameters", "psf_alpha2"]
        beta2 = block["psf_parameters", "psf_beta2"]
        alpha4 = block["psf_parameters", "psf_alpha4"]
        beta4 = block["psf_parameters", "psf_beta4"]
        
        uncorrelated_p_arr = np.array([alpha2, beta2, alpha4, beta4])
        p_arr = (np.linalg.inv(self.psf_transform_matrix)@uncorrelated_p_arr.T).T + self.psf_means
        delta_cl = np.zeros(len(self.psf_template[0][0]))
        for i in range(4):
            for j in range(4):
                delta_cl += p_arr[i]*p_arr[j]*self.psf_template[i][j]

        return delta_cl


    def extract_spectrum_prediction(self, block, data_type):
        """
        This is based on the original SACC likelihood, but adds the PSF subtraction
        """
        s = self.sacc_data

        # TODO: support 3x2pt here
        if data_type == "cl_ee":
            section = "shear_cl"
            is_auto = True
        elif data_type == "cl_bb":
            section = "shear_cl_bb"
            is_auto = True
        elif (data_type == "cl_eb") or (data_type == "cl_be"):
            section = "shear_cl_eb"
            is_auto = False
        else:
            raise ValueError(f"SACC likelihood does not yet understand data type {data_type}")

        ell_theory = block[section, "ell"]

        # We build up these vectors from all the data points.
        # Only the theory vector is needed for the likelihood - the others
        # are for convenience, debugging, etc.
        theory_vector = []
        angle_vector = []
        bin1_vector = []
        bin2_vector = []

        psf = self.compute_psf_template(block)

        # Because we called to_canonical_order when we loaded the data,
        # we know that the data is grouped by data type, and then by tracers (tomo bins).
        # So that means we can do a data type at a time and then concatenate them, and
        # within this do a bin pair at a time, and concatenate them too.
        for b1, b2 in s.get_tracer_combinations(data_type):
            # Here we assume that the bin names are formatted such that
            # they always end with _1, _2, etc. That isn't always true in
            # sacc, but is somewhat baked into cosmosis in other modules.
            # It would be nice to update that throughout, but that will
            # have to wait. Also, cosmosis bins start from 1 not 0.
            # We need to make sure that's fixed in the 
            i = int(b1.split("_")[-1]) + 1
            j = int(b2.split("_")[-1]) + 1


            try:
                cl_theory = block[section, f"bin_{i}_{j}"]
            except BlockError:
                if is_auto:
                    cl_theory = block[section, f"bin_{j}_{i}"]
                else:
                    raise

            # check that all the data points share the same window
            # object (window objects contain weights for a set of ell values,
            # as a matrix).
            window = None
            for d in s.get_data_points(data_type, (b1, b2)):
                w = d['window']
                if (window is not None) and (w is not window):
                    raise ValueError("Sacc likelihood currently assumes data types share a window object")
                window = w

            if window is None:
                raise ValueError("HSC likelihood assumes a window")
            # We need to interpolate between the sample ell values
            # onto all the ell values required by the weight function
            # This will give zero outside the range where we have
            # calculated the theory
            cl_theory_spline = sacc_like.SpectrumInterp(ell_theory, cl_theory)
            ell_window = window.values
            cl_theory_interpolated = cl_theory_spline(ell_window)

            for d in s.get_data_points(data_type, (b1, b2)):
                index = d['window_ind']
                ell_nominal = d['ell']
                weight = window.weight[:, index]

                # The weight away should hopefully sum to 1 anyway but we should
                # probably not rely on that always being true.
                binned_cl_theory = (weight @ cl_theory_interpolated) / weight.sum() + psf[index]

                theory_vector.append(binned_cl_theory)
                angle_vector.append(ell_nominal)
                bin1_vector.append(i - 1)
                bin2_vector.append(j - 1)

        # Return the whole collection as a single array
        theory_vector = np.array(theory_vector)

        # For convenience we also save the angle vector (ell or theta)
        # and bin indices
        angle_vector = np.array(angle_vector)
        bin1_vector = np.array(bin1_vector, dtype=int)
        bin2_vector = np.array(bin2_vector, dtype=int)

        return theory_vector, angle_vector, bin1_vector, bin2_vector


setup, execute, cleanup = HSCLikelihood.build_module()
