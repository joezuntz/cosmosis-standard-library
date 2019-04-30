from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
    return 0

def execute(block, config):
    A = block[names.halo_model_parameters,'A']
    eta0 = 1.03-0.11*A
    block[names.halo_model_parameters,'eta_0'] = eta0
    return 0

def cleanup(config):
    pass
