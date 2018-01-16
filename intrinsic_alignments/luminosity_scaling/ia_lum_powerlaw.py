# Based on the luminosity scaling in Joachimi et al (2011)
# Requires an interpolation table as the redshift scaling depends
# on the limiting magnitude and the galaxy sample in question

from cosmosis.datablock import names, option_section
import numpy as np


def setup(options):
    survey = block[options, 'survey']
    data = block[options, 'data']
    option = block[options, 'mode']

    if option not in ["luminosity_function", "interpolation_table"]:
        raise ValueError(
            "Please choose a method from the following: 'luminosity_function', ' interpolation_table'")

    if (option == "interpolation_table"):
        data = block[options, 'data']
    else:
        data = None

    config = survey, data
    return config


def execute(block, config):

    survey, filename = config
    ia_section = names.intrinsic_alignment_parameters

    # read in power spectra
    z, k, p_ii = block.get_grid(ia_section, "z", "k_h", "P_II")
    z, k, p_gi = block.get_grid(ia_section, "z", "k_h", "P_GI")

    beta = block[ia_section, 'beta']
    table = np.loadtxt(filename)

    L = table.T[2]
    z = table.T[0]
    rlim = table.T[1]

    grid = np.meshgrid(z, rlim)

    # Construct and apply luminosity scaling
    z_scaling = (1 + L)**beta
    p_ii *= z_scaling
    p_gi *= z_scaling

    # Save grid
    block.replace_grid(ia_section, "z", z, "k_h", k, "P_GI", p_gi)
    block.replace_grid(ia_section, "z", z, "k_h", k, "P_II", p_ii)
    return 0


def cleanup(config):
    pass
