"""

Adapted from MontePython 2.1, which is in turn adapted from 
arxiv: 1401.0358v2

"""
from __future__ import print_function
from cosmosis.datablock import names, option_section
import os
import numpy as np

dirname = os.path.split(__file__)[0]
default_data_file = os.path.join(dirname, "wigglez_bao.txt")
default_weight_file = os.path.join(dirname, "wigglez_bao_weight.txt")
default_rs_fiducial = 148.6


def setup(options):
    # Load the data from the default location unless otherwise specified
    data_file = options.get_string(
        option_section, "data_file", default_data_file)
    weight_file = options.get_string(
        option_section, "weight_file", default_weight_file)
    rs_fiducial = options.get_double(
        option_section, "rs_fiducial", default_rs_fiducial)
    verbose = options.get_bool(option_section, "verbose", False)

    # Load files
    z, Dv = np.loadtxt(data_file).T
    weight_matrix = np.loadtxt(weight_file)

    # Return data for later
    return (z, Dv, weight_matrix, rs_fiducial, verbose)


def execute(block, config):
    # Unpack data loaded in depending on options
    z_data, dv_data, weight_matrix, rs_fiducial, verbose = config

    # load distance relations and R_S
    z = block[names.distances, "z"]
    da = block[names.distances, "D_A"]  # in Mpc
    H = block[names.distances, "H"]  # in c/Mpc
    rs = block[names.distances, "RS_ZDRAG"]

    # Compute the derived D_V distance
    dr = z / H  # in Mpc/c
    dv = (da**2 * (1 + z)**2 * dr)**(1. / 3.)

    # Interpolate into the theory at the
    # observed redshifts
    dv_predicted = np.interp(z_data, z, dv)

    # If required, print out some info
    if verbose:
        print("dv_predicted = ",  dv_predicted)
        print("dv_data = ", dv_data)
        print("rs = ", rs)
        print("rs_fiducial = ", rs_fiducial)

    # Get the Gaussian likelihood, including the scaling wrt the fiducial
    delta = dv_data - dv_predicted / rs * rs_fiducial
    like = -0.5 * np.einsum('i,ij,j', delta, weight_matrix, delta)

    # Save the result
    block[names.likelihoods, "wigglez_bao_like"] = like

    # Signal success
    return 0


def cleanup(config):
    pass
