from cosmosis.datablock import names, option_section
from luminosity_function import jb_calculate_alpha
import luminosity_function as luminosity
import numpy as np


def setup(options):
    # method = options[option_section, "method"].lower()
    # For the moment we only have one method. It may or may not be useful to try implementing more here later.
    # if method not in ["jb"]:
    #	raise ValueError('The method in the luminosity function module must'
    #		'be either "JB" (Joachimi and Bridle 2010) or')
    method = "jb"

    try:
        mag_lim = options[option_section, "magnitude_limit"]
    except:
        mag_lim = 24
    try:
        binned_alpha = options[option_section, "binned_alpha"]
    except:
        binned_alpha = True
    try:
        sample = options[option_section, "sample"]
    except:
        sample = 'wl_number_density'

    config = {'method': method, 'mag_lim': mag_lim,
              'use_binned_alpha': binned_alpha,
              'sample':sample
              }
    return config


def execute(block, config):
    ia = names.intrinsic_alignment_parameters
    cos = names.cosmological_parameters
    lum = names.galaxy_luminosity_function

    method = config['method']
    mag_lim = config['mag_lim']
    sample = config['sample']

    # Obtain the fine grid points and limits for the redshift from the datablock
    # If the information isn't there already, set these parameters to sensible values
    try:
        Nz = block[sample, 'nz']
        nzbin = block[sample, 'nbin']
        zmax = block[sample, 'edge_%d' % (nzbin + 1)]
    except:
        Nz = 500
        zmax = 3.0

    coeff_a = luminosity.initialise_jb_coefficients(mag_lim)
    alpha, z = jb_calculate_alpha(coeff_a, zmax, Nz)

    # Write the fine grid alpha to the datablock
    block.put_double_array_1d(lum, 'alpha', alpha)
    block.put_double_array_1d(lum, 'z', z)

    # If required then interpolate alpha(z,rlim) to the mean values in each of the specified redshift bins
    use_binned_alpha = config['use_binned_alpha']
    if use_binned_alpha:
        alpha_bin, z_bar = luminosity.get_binned_alpha(block, alpha, z, sample)

        # Then write these to the datablock
        block.put_double_array_1d(lum, 'alpha_binned', alpha_bin)
        block.put_double_array_1d(lum, 'z_binned', z_bar)

    return 0


def cleanup(config):
    pass
