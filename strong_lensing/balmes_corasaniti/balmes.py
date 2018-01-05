from cosmosis.datablock import names, option_section
import numpy as np
import os

DIRNAME = os.path.split(__file__)[0]
DEFAULT_DATA_FILE = os.path.join(DIRNAME, "balmes.txt")


def setup(options):
    data_file = options.get_string(
        option_section, "data_file", default=DEFAULT_DATA_FILE)
    H0, P = np.loadtxt(data_file).T
    return [H0, P]


def execute(block, config):
    H0_data, P_data = config
    H0_theory = block[names.cosmological_parameters, "h0"] * 100
    if (H0_theory < H0_data[0]) or (H0_theory > H0_data[-1]):
        logP = -np.inf
    else:
        P = np.interp(H0_theory, H0_data, P_data)
        logP = np.log(P)
    block[names.likelihoods, "balmes_sl_like"] = logP
    return 0
