"""
An interface to the baryonic code.

"""
from __future__ import print_function
from builtins import range
from cosmosis.datablock import option_section, names as section_names
import baryonic
import numpy as np
import os

baryonPowerModulator = None
ALREADY_DONE_KEY = "BARYONICEFFECTSDONE"


def setup(options):
    section = option_section
    mode = options.get_string(option_section, "mode", "ratio")

    try:
        ratiotable = options[option_section, 'ratiotable']
    except KeyError:
        raise ValueError(
            "Please set the parameter powtable in the owls section of the ini file, pointing to a ratio file to use")
        #print("Selected ratio-mode baryons.  Loading table from ", ratiotable)

    baryonPowerModulator = baryonic.RatioDatPowerModulator(ratiotable)

    return baryonPowerModulator


"""
def load_parameters(package, modulator):
	section = pydesglue.section_names.baryon_parameters
	if isinstance(modulator, owls.ChebyshevBaryonPowerModulator):
		rk = [package.get_param(section, "OWLS_K%d"%i) for i in xrange(1,modulator.nterm+1)]
		rz = [package.get_param(section, "OWLS_Z%d"%i) for i in xrange(1,modulator.nterm+1)]
		r = np.concatenate([rk, rz])
	elif isinstance(modulator, owls.FixedBaryonPowerModulator):
		r=None
	elif isinstance(modulator, owls.BaryonPowerModulator):
		r = package.get_param(section, "OWLS_R")
	elif isinstance(modulator, owls.ScaledBaryonPowerModulator):
		r = package.get_param(section, "OWLS_R")	
	else:
		raise ValueError("Internal error - modulator type not understood in OWLS")
	return r
"""


def load_parameters(block, modulator):
    r=None
    return r


def execute(block, config):
    # JAZ some day we will switch to get this from an input,
    # so that the modulator only happens once
    modulator = config
    r = load_parameters(block, modulator)
    section = section_names.matter_power_nl
    # load other current values
    k, z, P = (block.get_grid(section, "k_h", "z", "P_k"))
    print(P.shape)
    P_mod = modulator.modulate(k, z, P, r)
    block.replace_grid(section, "k_h", k, "z", z, "P_k", P_mod)
    return 0
