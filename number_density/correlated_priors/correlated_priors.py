import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    input_parameters = options.get_string(option_section, "uncorrelated_parameters").split()
    input_parameters = [p.strip().split("/") for p in input_parameters]
    output_parameters = options.get_string(option_section, "output_parameters").split()
    output_parameters = [p.strip().split("/") for p in output_parameters]
    covariance_file = options[option_section, "covariance"]

    cov = np.loadtxt(covariance_file)
    L = np.linalg.cholesky(cov) 
    return input_parameters, output_parameters, L

def execute(block, config):
    input_parameters, output_parameters, L = config
    
    p = []
    for section, name in input_parameters:
        p.append(block[section, name])
    p = L @ np.array(p)

    for i, (section, name) in enumerate(output_parameters):
        block[section, name] = p[i]

    return 0

def clean(config):
    pass
