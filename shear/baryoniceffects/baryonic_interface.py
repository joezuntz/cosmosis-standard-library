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
    if mode == "fixed":
        try:
            powtable = options[option_section, 'powtable']
            dmtable = options[option_section, 'dmtable']
        except KeyError:
            raise ValueError(
                "Please set the parameter powtable in the baryonic section of the ini file, pointing to an powtable file to use")
        print("Selected fixed-mode baryons.  Loading table from ", powtable,dmtable)
        baryonPowerModulator = baryonic.FixedBaryonPowerModulator(powtable,dmtable)
        
    elif mode == "ratio":
        try:
            ratiotable = options[option_section, 'ratiotable']
        except KeyError:
            raise ValueError(
                "Please set the parameter powtable in the owls section of the ini file, pointing to a ratio file to use")
        print("Selected ratio-mode baryons.  Loading table from ", ratiotable)
        ftype = ratiotable.split('.')[-1]
        if ftype == 'fits':
            baryonPowerModulator = baryonic.RatioTablePowerModulator(ratiotable)
        elif ftype == 'dat':
            baryonPowerModulator = baryonic.RatioDatPowerModulator(ratiotable)
        else:
            print(ratiotable, "should have extension fits or dat.")

    else:
        raise ValueError(
            "'mode' must be either 'fixed' (no params; fixed OWLS sim), or 'ratio'.")
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
    section = "baryonic"
    if isinstance(modulator, baryonic.ChebyshevBaryonPowerModulator):
        rk = [block[section, "K%d" % i]
              for i in range(1, modulator.nterm + 1)]
        rz = [block[section, "Z%d" % i]
              for i in range(1, modulator.nterm + 1)]
        r = np.concatenate([rk, rz])
    elif isinstance(modulator, baryonic.FixedBaryonPowerModulator):
        r = None
    elif isinstance(modulator, baryonic.RatioTablePowerModulator):
        r = None
    elif isinstance(modulator, baryonic.RatioDatPowerModulator):
        r = None
    elif isinstance(modulator, baryonic.BaryonPowerModulator):
        r = block[section, "R"]
    elif isinstance(modulator, baryonic.ScaledBaryonPowerModulator):
        r = block[section, "R"]
    else:
        raise ValueError(
            "Internal error - modulator type not understood in OWLS")
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
