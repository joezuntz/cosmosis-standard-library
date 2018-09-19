"""
This module saves cosmosis output to a 2-pt file after interpolating it 
to specified ell values. It is useful for making simulations, but
is not yet fully supported so please use it with care and check 
the results carefully.

"""
from __future__ import print_function


from builtins import range
from builtins import object
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
import twopoint
from twopoint_cosmosis import type_table
import gaussian_covariance
twopt_like = __import__('2pt_like')  # don't start .py files with a number!
SpectrumInterp = twopt_like.SpectrumInterp


def setup(options):
    # ell range of output - all assumed the same, with log-spacing and
    # no window functions
    real_space = options.get_bool(option_section, "real_space", False)
    make_covariance = options.get_bool(option_section, "make_covariance", True)
    

    config = {"real_space": real_space, "make_covariance": make_covariance}

    if real_space:
        if options.has_value(option_section, "theta"):
            theta = options.get_double_array_1d(option_section, "theta")
        else:
            theta_min = options.get_double(option_section, "theta_min")
            theta_max = options.get_double(option_section, "theta_max")
            n_theta = options.get_int(option_section, "n_theta")
            theta = np.logspace(np.log10(theta_min), np.log10(
                theta_max), n_theta + 1, endpoint=True)
            a = theta[:-1]
            b = theta[1:]
            theta = 2. / 3. * (b**3 - a**3) / (b**2 - a**2)
        config['theta'] = theta
        print("Saving at these theta values:")
        print(theta)
    else:
        if options.has_value(option_section, "ell"):
            ell = options.get_double_array_1d(option_section, "ell")
        else:
            ell_min = options.get_double(option_section, "ell_min")
            ell_max = options.get_double(option_section, "ell_max")
            n_ell = options.get_int(option_section, "n_ell")
            ell = np.logspace(np.log10(ell_min), np.log10(ell_max), n_ell)
        config['ell'] = ell
        print("Saving at these ell values:")
        print(ell)

    def get_arr(x):
        a = options[option_section, x]
        if not isinstance(a, np.ndarray):
            a = np.array([a])
        return a

    if make_covariance:
        config['number_density_shear_bin'] = get_arr(
            "number_density_shear_bin")
        config['number_density_lss_bin'] = get_arr("number_density_lss_bin")
        config['sigma_e_bin'] = get_arr("sigma_e_bin")
        config['survey_area'] = options[option_section,
                                        "survey_area"] * (np.pi * np.pi) / (180 * 180)

    # names n(z) sections in the datablock, to be saved with the same name
    # to the FITS file output
    config['shear_nz'] = options.get_string(option_section, "shear_nz_name").upper()
    config['position_nz'] = options.get_string(option_section, "position_nz_name").upper()
    # Whether to add nz_ to the start of these
    config['prefix_nz_section'] = options.get_bool(option_section, "prefix_nz_section", True)

    # name of the output file and whether to overwrite it.
    config['filename'] = options.get_string(option_section, "filename")
    config['clobber'] = options.get_bool(option_section, "clobber", False)

    scale_cuts = {}
    bin_cuts = []
    range_token = "angle_range_"
    cut_token = "cut_"
    for _, key in options.keys(option_section):
        if key.startswith(range_token):
            bits = key[len(range_token):].split("_")
            name = "_".join(bits[:-2])
            bin1 = int(bits[-2])
            bin2 = int(bits[-1])
            scale_cuts[(name, bin1, bin2)] = options.get_double_array_1d(
                option_section, key)
        elif key.startswith(cut_token):
            name = key[len(cut_token):]
            cuts = options[option_section, key].split()
            cuts = [eval(cut) for cut in cuts]
            for b1, b2 in cuts:
                bin_cuts.append((name, b1, b2))
    config['bin_cuts'] = bin_cuts
    config['scale_cuts'] = scale_cuts

    return config


def spectrum_measurement_from_block(block, section_name, output_name, types, kernels, angle_sample, real_space):

    # The dictionary type_table stores the codes used in the FITS files for
    # the types of spectrum
    type_codes = (types[0].name, types[1].name)
    _, _, bin_format = type_table[type_codes]

    # for cross correlations we must save bin_ji as well as bin_ij.
    # but not for auto-correlations. Also the numbers of bins can be different
    is_auto = (types[0] == types[1])
    if block.has_value(section_name, "nbin_a"):
        nbin_a = block[section_name, "nbin_a"]
        nbin_b = block[section_name, "nbin_b"]
    else:
        nbin_a = block[section_name, "nbin"]
        nbin_b = block[section_name, "nbin"]

    # This is the ell values that have been calculated by cosmosis, not to
    # be confused with the ell values at which we want to save the results
    #(which is ell_sample)
    if real_space:
        # This is in radians
        theory_angle = block[section_name, "theta"]
        angle_sample_radians = np.radians(angle_sample / 60.)
    else:
        theory_angle = block[section_name, "ell"]
        angle_sample_radians = angle_sample

    # This is the length of the sample values
    n_angle_sample = len(angle_sample)

    # The fits format stores all the measurements
    # as one long vector.  So we build that up here from the various
    # bins that we will load in.  These are the different columns
    bin1 = []
    bin2 = []
    value = []
    angular_bin = []
    angle = []

    # Bin pairs. Varies depending on auto-correlation
    for i in range(nbin_a):
        if is_auto:
            jstart = i
        else:
            jstart = 0
        for j in range(jstart, nbin_b):
            if is_auto:
                # Load and interpolate from the block
                bin_name = bin_format.format(i + 1, j + 1)
                if block.has_value(section_name, bin_name):
                    cl = block[section_name, bin_format.format(i + 1, j + 1)]
                else:
                    cl = block[section_name, bin_format.format(j + 1, i + 1)]
            else:
                cl = block[section_name, bin_format.format(i + 1, j + 1)]
            # Convert arcmin to radians for the interpolation
            cl_interp = SpectrumInterp(theory_angle, cl)
            cl_sample = cl_interp(angle_sample_radians)
            #cl_sample = interp1d(theory_angle, cl)(angle_sample_radians)
            # Build up on the various vectors that we need
            bin1.append(np.repeat(i + 1, n_angle_sample))
            bin2.append(np.repeat(j + 1, n_angle_sample))
            value.append(cl_sample)
            angular_bin.append(np.arange(n_angle_sample))
            angle.append(angle_sample)

    # Convert all the lists of vectors into long single vectors
    bin1 = np.concatenate(bin1)
    bin2 = np.concatenate(bin2)
    angular_bin = np.concatenate(angular_bin)
    value = np.concatenate(value)
    angle = np.concatenate(angle)
    bins = (bin1, bin2)

    # At the moment we only support this window function
    windows = "SAMPLE"

    if real_space:
        angle_unit = "arcmin"
    else:
        angle_unit = None
    # Build the output object type reqired.
    s = twopoint.SpectrumMeasurement(output_name, bins, types, kernels, windows,
                                     angular_bin, value, angle=angle, angle_unit=angle_unit)
    return s




def nz_from_block(block, nz_name, prefix_nz_section):
    print()
    print()
    print("*************************************************************************************")
    print("Saving n(z) from the block to file.")
    print("A quick warning - we are assuming things about the n(z) that may not be quite right.")
    print("Converting from splines to histograms.")
    print("To properly fix this I will have to do a bit more work.")
    print()
    print("*************************************************************************************")
    print("ANOTHER WARNING: this has recently changed to hopefully be slightly better")
    print("It may be different to your previous results.")
    print("*************************************************************************************")
    print()
    print()
    
    if prefix_nz_section:
        section_name = "nz_"+nz_name
    else:
        section_name = nz_name

    z = block[section_name, "z"]
    dz = 0.5 * (z[10] - z[9])
    zlow = z - dz
    zhigh = z + dz
    if zlow[0] < 0:
        zlow = zlow[1:]
        z = z[1:]
        zhigh = zhigh[1:]
        cut = True
    else:
        cut = False
    assert zlow[0] > 0
    nbin = block[section_name, "nbin"]
    nzs = []
    for i in range(nbin):
        nz = block[section_name, "bin_{}".format(i + 1)]
        if cut:
            nz = nz[1:]
        nzs.append(nz)

    return twopoint.NumberDensity(nz_name, zlow, z, zhigh, nzs)


def convert_nz_steradian(n):
    return n * (41253.0 * 60. * 60.) / (4 * np.pi)


class ObservedClGetter(object):
    def __init__(self, number_density_shear_bin, number_density_lss_bin, sigma_e_bin):
        self.number_density_shear_bin = convert_nz_steradian(
            number_density_shear_bin)
        self.number_density_lss_bin = convert_nz_steradian(
            number_density_lss_bin)
        self.sigma_e_bin = sigma_e_bin
        self.splines = {}

    def lookup(self, block, A, B, i, j, ell):
        """
        This is a helper function for the compute_gaussian_covariance code.
        It looks up the theory value of C^{ij}_{AB}(ell) in the block and adds noise
        """
        # We have already saved splines into the theory space earlier
        # when constructing the theory vector.
        # So now we just need to look those up again, using the same
        # code we use in the twopoint library.
        section, ell_name, value_name = type_table[A, B]
        assert ell_name == "ell", "Gaussian covariances are currently only written for C_ell, not other 2pt functions"

        name_ij = value_name.format(i, j)
        section_name_ij = '{}_{}'.format(section, name_ij)

        if name_ij in self.splines:
            spline = self.splines[section_name_ij]
        else:
            spline = self.make_spline(block, A, B, i, j, ell)
            self.splines[section_name_ij] = spline

        obs_c_ell = spline(ell)

        # For shear-shear the noise component is sigma^2 / number_density_bin
        # and for position-position it is just 1/number_density_bin
        if (A == B) and (A == twopoint.Types.galaxy_shear_emode_fourier.name) and (i == j):
            noise = self.sigma_e_bin[i - 1]**2 / \
                self.number_density_shear_bin[i - 1]
            obs_c_ell += noise
        if (A == B) and (A == twopoint.Types.galaxy_position_fourier.name) and (i == j):
            noise = 1.0 / self.number_density_lss_bin[i - 1]
            obs_c_ell += noise

        return obs_c_ell

    def make_spline(self, block, A, B, i, j, ell):
        section, ell_name, value_name = type_table[A, B]
        assert ell_name == "ell", "Gaussian covariances are currently only written for C_ell, not other 2pt functions"

        # We extract relevant bits from the block and spline them
        # for output
        name_ij = value_name.format(i, j)
        name_ji = value_name.format(j, i)

        angle = block[section, ell_name]

        if block.has_value(section, name_ij):
            theory = block[section, name_ij]
        elif block.has_value(section, name_ji) and A == B:
            theory = block[section, name_ji]
        else:
            raise ValueError("Could not find theory prediction {} in section {}".format(name_ij, section))

        spline = interp1d(angle, theory)
        return spline


def covmat_from_block(block, spectra, sky_area, number_density_shear_bin, number_density_lss_bin, sigma_e_bin):
    getter = ObservedClGetter(number_density_shear_bin,
                              number_density_lss_bin, sigma_e_bin)
    C = []
    names = []
    starts = []
    lengths = []

    # s and t index the spectra that we have. e.g. s or t=1 might be the full set of
    # shear-shear measuremnts
    x = 0
    for s, AB in enumerate(spectra[:]):
        M = []
        starts.append(x)
        L = len(AB)
        lengths.append(L)
        x += L

        names.append(AB.name)
        for t, CD in enumerate(spectra[:]):
            print("Looking at covariance between {} and {} (s={}, t={})".format(AB.name, CD.name, s, t))
            # We only calculate the upper triangular.
            # Get the lower triangular here. We have to
            # transpose it compared to the upper one.
            if s > t:
                MI = C[t][s].T
            else:
                MI = gaussian_covariance.compute_gaussian_covariance(sky_area,
                                                                     getter.lookup, block, AB, CD)
            M.append(MI)
        C.append(M)
    C = np.vstack([np.hstack(CI) for CI in C])
    info = twopoint.CovarianceMatrixInfo("COVMAT", names, lengths, C)
    return info


def execute(block, config):
    real_space = config['real_space']
    fourier_space = not real_space
    if real_space:
        angle_sample = config['theta']
    else:
        angle_sample = config['ell']

    filename = config['filename']
    shear_nz = config['shear_nz']
    position_nz = config['position_nz']
    clobber = config['clobber']

    make_covariance = config['make_covariance']
    if make_covariance:
        number_density_shear_bin = config['number_density_shear_bin']
        number_density_lss_bin = config['number_density_lss_bin']
        sigma_e_bin = config['sigma_e_bin']
        survey_area = config['survey_area']

    print("Saving two-point data to {}".format(filename))

    spectra = []

    need_source_nz = False
    need_lens_nz = False
    if fourier_space and block.has_section(names.shear_cl):
        name = "shear_cl"
        types = (twopoint.Types.galaxy_shear_emode_fourier,
                 twopoint.Types.galaxy_shear_emode_fourier)
        kernels = (shear_nz, shear_nz)
        s = spectrum_measurement_from_block(
            block, name, name, types, kernels, angle_sample, real_space)
        print(" - saving shear_cl")
        need_source_nz = True
        spectra.append(s)

    if fourier_space and block.has_section("galaxy_shear_cl"):
        name = "galaxy_shear_cl"
        types = (twopoint.Types.galaxy_position_fourier,
                 twopoint.Types.galaxy_shear_emode_fourier,)
        kernels = (position_nz, shear_nz)
        s = spectrum_measurement_from_block(
            block, name, name, types, kernels, angle_sample, real_space)
        print(" - saving galaxy_shear_cl")
        need_source_nz = True
        need_lens_nz = True
        spectra.append(s)

    if fourier_space and block.has_section("galaxy_cl"):
        name = "galaxy_cl"
        types = (twopoint.Types.galaxy_position_fourier,
                 twopoint.Types.galaxy_position_fourier)
        kernels = (position_nz, position_nz)
        s = spectrum_measurement_from_block(
            block, name, name, types, kernels, angle_sample, real_space)
        print(" - saving galaxy_cl")
        need_lens_nz = True
        spectra.append(s)

    if real_space and block.has_section("shear_xi"):
        types = (twopoint.Types.galaxy_shear_plus_real,
                 twopoint.Types.galaxy_shear_plus_real)
        kernels = (shear_nz, shear_nz)
        s = spectrum_measurement_from_block(
            block, "shear_xi", "xip", types, kernels, angle_sample, real_space)
        print(" - saving xi_plus")
        need_source_nz = True
        spectra.append(s)

        name = "shear_xi"
        types = (twopoint.Types.galaxy_shear_minus_real,
                 twopoint.Types.galaxy_shear_minus_real)
        kernels = (shear_nz, shear_nz)
        s = spectrum_measurement_from_block(
            block, "shear_xi", "xim", types, kernels, angle_sample, real_space)
        need_source_nz = True
        print(" - saving xi_minus")
        spectra.append(s)

    if real_space and block.has_section("galaxy_shear_xi"):
        types = (twopoint.Types.galaxy_position_real,
                 twopoint.Types.galaxy_shear_plus_real)
        kernels = (position_nz, shear_nz)
        s = spectrum_measurement_from_block(
            block, "galaxy_shear_xi", "gammat", types, kernels, angle_sample, real_space)
        print(" - saving gammat")
        need_source_nz = True
        need_lens_nz = True
        spectra.append(s)

    if real_space and block.has_section("galaxy_xi"):
        types = (twopoint.Types.galaxy_position_real,
                 twopoint.Types.galaxy_position_real)
        kernels = (position_nz, position_nz)
        s = spectrum_measurement_from_block(
            block, "galaxy_xi", "wtheta", types, kernels, angle_sample, real_space)
        print(" - saving wtheta")
        need_lens_nz = True
        spectra.append(s)

    if make_covariance:
        covmat_info = covmat_from_block(
            block, spectra, survey_area, number_density_shear_bin, number_density_lss_bin, sigma_e_bin)
    else:
        covmat_info = None

    if not spectra:
        raise ValueError("Sorry - I couldn't find any spectra to save.")

    kernels = []
    if need_source_nz:

        kernels.append(nz_from_block(block, shear_nz, config['prefix_nz_section']))
    if need_lens_nz and (shear_nz != position_nz):
        kernels.append(nz_from_block(block, position_nz, config['prefix_nz_section']))

    windows = []

    data = twopoint.TwoPointFile(spectra, kernels, windows, covmat_info)

    # Apply cuts
    scale_cuts = config['scale_cuts']
    bin_cuts = config['bin_cuts']
    if scale_cuts or bin_cuts:
        data.mask_scales(scale_cuts, bin_cuts)

    data.to_fits(filename, clobber=clobber)

    return 0
