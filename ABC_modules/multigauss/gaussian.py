from builtins import str
from builtins import range
import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import os


cosmo = section_names.cosmological_parameters

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]
COV_file = os.path.join(ROOT_dir, 'abc_multigauss_cov.txt')


def setup(options):
    section = option_section
    ngauss = options.get_int(section, "ngauss", default=4)
    covfile = options.get_string(section, "cov_file", default=COV_file)
    ndata = options.get_int(section, "ndata", default=5000)
    sigma = np.loadtxt(covfile)

    return (sigma, ndata)


def execute(block, config):

    # Simple ABC model which outputs ndata simulated data points from an ngauss multigaussian

    sigma, ndata = config

    params = np.zeros(len(sigma[0]))
    for i in range(len(params)):
        params[i] = block[cosmo, 'mu' + str(i)]

    model = np.random.multivariate_normal(params, sigma, ndata)
    block[section_names.data_vector, 'abc_multigauss_simulation'] = model

    return 0


def cleanup(config):
    return 0
