from cosmosis.datablock import names, option_section
import os
import numpy as np


dirname = os.path.split(__file__)[0]


def setup(options):
    mass = options[option_section, "mass"]
    if mass not in ['low', 'high']:
        raise ValueError("Please choose low or high mass for Jullo likelihood")
    data = np.loadtxt(os.path.join(dirname, "jullo_data.txt")).T
    if mass == "low":
        data = data[0], data[1], data[2]
    else:
        data = data[0], data[3], data[4]
    return data, mass


def execute(block, config):
    data, mass = config
    z_obs, b_obs, sigma_obs = data
    # Load the bias from the block
    try:
        z = block[names.bias_field, "z"]
    except:
        raise ValueError(
            "The Jullo data requires bias as a function of redshift")
    b = block[names.bias_field, "b"]

    # Just use smallest k value for 2D
    if b.ndim == 2:
        k, z, b = block.get_grid(names.bias_field, "k_h", "z", "b")
        b = b[0]

    # interpolate into it at the data z
    b_theory = np.interp(z_obs, z, b)

    # get likelihood
    chi2 = ((b_theory - b_obs)**2 / sigma_obs**2).sum()
    like = -0.5 * chi2
    block[names.likelihoods, "jullo_like"] = like

    return 0


def cleanup(config):
    pass
