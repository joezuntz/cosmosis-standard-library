"""
An interface to the OWLS code.

"""
from __future__ import print_function
from builtins import range
from cosmosis.datablock import option_section, names as section_names
import owls
import numpy as np
import os

baryonPowerModulator = None
ALREADY_DONE_KEY = "OWLSDONE"


def setup(options):
    section = option_section
    mode = options.get_string(option_section, "mode", "single")
    if mode == "single":
        baryonPowerModulator = owls.BaryonPowerModulator()
        print("Using single-parameter OWLS mode")
    elif mode == "chebyshev":
        try:
            nterm = options['nterm']
        except KeyError:
            raise ValueError(
                "Please set the parameter nterm in the owls section of the ini file to use Chebyshev OWLS")
        baryonPowerModulator = owls.ChebyshevBaryonPowerModulator(
            nterm=nterm, extremes=1.0)
        print("Using checbyshev mode with %d terms" % nterm)
    elif mode == "fixed":
        try:
            powtable = options[option_section, 'powtable']
        except KeyError:
            raise ValueError(
                "Please set the parameter powtable in the owls section of the ini file, pointing to an OWLS powtable file to use")
        print("Selected fixed-mode baryons.  Loading table from ", powtable)
        baryonPowerModulator = owls.FixedBaryonPowerModulator(powtable)
    elif mode == "scaled":
        try:
            powtable = options[option_section, 'powtable']
        except KeyError:
            raise ValueError(
                "Please set the parameter powtable in the owls section of the ini file, pointing to an OWLS powtable file to use")
        print("Selected scaled-mode baryons.  Loading table from ", powtable)
        baryonPowerModulator = owls.ScaledBaryonPowerModulator(powtable)

    else:
        raise ValueError(
            "'mode' must be either 'single' (one parameter, OWLS_R), 'chebyshev' (nterm parameters OWLS_R1,...), 'fixed' (no params; fixed OWLS sim), or 'scaled' (one parameter, OWLS_R) ")
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
    section = "OWLS"
    if isinstance(modulator, owls.ChebyshevBaryonPowerModulator):
        rk = [block[section, "K%d" % i]
              for i in range(1, modulator.nterm + 1)]
        rz = [block[section, "Z%d" % i]
              for i in range(1, modulator.nterm + 1)]
        r = np.concatenate([rk, rz])
    elif isinstance(modulator, owls.FixedBaryonPowerModulator):
        r = None
    elif isinstance(modulator, owls.BaryonPowerModulator):
        r = block[section, "R"]
    elif isinstance(modulator, owls.ScaledBaryonPowerModulator):
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
