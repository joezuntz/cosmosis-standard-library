"""
An interface to the baryonic code.

"""
from __future__ import print_function
from cosmosis.datablock import option_section, names as section_names
import baryonic
import numpy as np
import os



def setup(options):
    section = option_section

    try:
        ratio_table = options[option_section, 'ratio_table']
    except KeyError:
        raise ValueError("Please set the parameter powtable in the owls section "
            "of the ini file, pointing to a ratio file to use")

    # This objejct does the actual work
    baryonPowerModulator = baryonic.RatioDatPowerModulator(ratio_table)

    return baryonPowerModulator


def execute(block, config):
    modulator = config

    # Load the current unmodulated matter power spectrum
    section = section_names.matter_power_nl
    k, z, P = (block.get_grid(section, "k_h", "z", "P_k"))
    
    # Apply the scaling as a function of k and z
    P_mod = modulator.modulate(k, z, P)

    # Update the stored values
    block.replace_grid(section, "k_h", k, "z", z, "P_k", P_mod)

    # All is good - return
    return 0
