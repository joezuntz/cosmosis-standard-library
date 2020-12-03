from __future__ import print_function
from cosmosis.datablock import option_section, names
import consistency
import numpy as np


def setup(options):
    verbose = options.get_bool(option_section, "verbose", default=False)
    cosmomc_theta = options.get_bool(option_section, "cosmomc_theta", default=False)
    relations_file = options.get_string(
        option_section, "relations_file", default="")
    cons = consistency.cosmology_consistency(verbose, relations_file, cosmomc_theta)
    return cons


def execute(block, config):
    cons = config
    cosmo = names.cosmological_parameters

    # Create dict of all parameters that we have already
    known_parameters = {}
    for param in cons.parameters:
        if '___' in param:
            section, lookup_param = param.split('___')
        else:
            section = cosmo
            lookup_param = param
        if block.has_value(section, lookup_param):
            known_parameters[param] = block[section, lookup_param]

    if cons.verbose:
        print("Consistency relation input parameters:", ', '.join(list(known_parameters.keys())))

    # Run the consistency checker/parameter filler-inner.
    try:
        filled_parameters = cons(known_parameters)
    # There are two possible simpler errors:
    # too many/inconsistent parameters, or not enough.
    except consistency.UnderSpecifiedModel as error:
        print("You did not set enough cosmological parameters")
        print("to be able to deduce the rest of them:")
        print(error)
        return 1
    except consistency.OverSpecifiedModel as error:
        print("You set inconsistent cosmological parameters:")
        print(error)
        return 2

    # Annoyingly this does not fit elsewhere
    if block.has_value(cosmo, "log1e10As"):
        block[cosmo, "A_s"] = np.exp(block[cosmo, 'log1e10As']) * 1.0e-10
    elif block.has_value(cosmo, "A_s_1e9"):
        block[cosmo, "A_s"] = block[cosmo, "A_s_1e9"] * 1e-9

    # Set or replace the new values
    for param, value in list(filled_parameters.items()):
        if '___' in param:
            section, param = param.split('___')
        else:
            section = cosmo
        block[section, param] = value

    return 0
