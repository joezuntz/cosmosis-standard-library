import numpy as np
from cosmosis.datablock import BlockError
import pathlib
import sys

# Get the SpectrumInterp class from the spec_tools module.
# Should really put this somewhere else!
twopoint_dir = pathlib.Path(__file__).parent.parent.parent.resolve() / "2pt"
print(twopoint_dir)
sys.path.append(str(twopoint_dir))
from spec_tools import SpectrumInterp


def extract_spectrum_prediction(sacc_data, block, data_type, section, **kwargs):


        ell_theory = block[section, "ell"]
        is_auto = block[section, "is_auto"]

        # We build up these vectors from all the data points.
        # Only the theory vector is needed for the likelihood - the others
        # are for convenience, debugging, etc.
        theory_vector = []
        angle_vector = []
        bin1_vector = []
        bin2_vector = []


        # Because we called to_canonical_order when we loaded the data,
        # we know that the data is grouped by data type, and then by tracers (tomo bins).
        # So that means we can do a data type at a time and then concatenate them, and
        # within this do a bin pair at a time, and concatenate them too.
        for b1, b2 in sacc_data.get_tracer_combinations(data_type):
            # Here we assume that the bin names are formatted such that
            # they always end with _1, _2, etc. That isn't always true in
            # sacc, but is somewhat baked into cosmosis in other modules.
            # It would be nice to update that throughout, but that will
            # have to wait. Also, cosmosis bins start from 1 not 0.
            # We need to make sure that's fixed in the 
            i = int(b1.split("_")[-1]) + 1
            j = int(b2.split("_")[-1]) + 1

            if data_type in kwargs.get("flip", False):
                i, j = j, i

            try:
                cl_theory = block[section, f"bin_{i}_{j}"]
            except BlockError:
                if is_auto:
                    cl_theory = block[section, f"bin_{j}_{i}"]
                else:
                    raise

            # check that all the data points share the same window
            # object (window objects contain weights for a set of ell values,
            # as a matrix), or that none have windows.
            window = None
            for d in sacc_data.get_data_points(data_type, (b1, b2)):
                w = d.get_tag('window')
                if (window is not None) and (w is not window):
                    raise ValueError("Sacc likelihood currently assumes data types share a window object")
                window = w

            # We need to interpolate between the sample ell values
            # onto all the ell values required by the weight function
            # This will give zero outside the range where we have
            # calculated the theory
            cl_theory_spline = SpectrumInterp(ell_theory, cl_theory)
            if window is not None:
                ell_window = window.values
                cl_theory_interpolated = cl_theory_spline(ell_window)

            for d in sacc_data.get_data_points(data_type, (b1, b2)):
                ell_nominal = d['ell']
                if window is None:
                    binned_cl_theory = cl_theory_spline(ell_nominal)
                else:
                    index = d['window_ind']
                    weight = window.weight[:, index]

                    # The weight away should hopefully sum to 1 anyway but we should
                    # probably not rely on that always being true.
                    binned_cl_theory = (weight @ cl_theory_interpolated) / weight.sum()

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

        metadata = {
            "angle": angle_vector,
            "bin1": bin1_vector,
            "bin2": bin2_vector
        }

        return theory_vector, metadata
