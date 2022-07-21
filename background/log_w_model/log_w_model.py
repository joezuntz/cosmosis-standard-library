"""
This module, which is mainly here as an example, implements
one w(z) module used in Tripathi, Sangwan, Jassal (2017), in which
w(z) = w_0 + w_1 log(1 + z)

It saves to the de_equation_of_state section, which can be read by camb.
"""
from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
    zmin = 0.0

    # These options are read from the first parameter file,
    # and are only read once at the beginning of the entire
    # analysis - they are fixed over the course of a chain.
    zmax = options.get_double(option_section, "zmax", 3.0)
    nz = options.get_int(option_section, "nz", 301)

    # anything we return here will be saved and be sent
    # to the execute function below as the second argument.
    z = np.linspace(zmin, zmax, nz)[::-1]
    a = 1 / ( 1 + z)
    log_1_plus_z = np.log(1 + z)
    return {"z": z, "a": a, "log_term": log_1_plus_z}


def execute(block, config):
    # The block argument is values that have been either supplied
    # by the sampler or computed by earlier stages of the pipeline.
    # The config argument is whatever the setup function above saved
    # earlier.
    z = config["z"]
    a = config["a"]
    log_1_plus_z = config["log_term"]


    # We load things made by earlier pipeline stages into the block like this.
    # The two keys are both strings (the names module just contains a collection
    # of standard category names, but we can use anything we want).
    w0 = block[names.cosmological_parameters, "w"]
    w1 = block[names.cosmological_parameters, "w1"]

    # The actual calculation
    w = w0 + w1 * log_1_plus_z

    # Now we save our results back to the data block for future pipeline stages
    # to use. CAMB is looking for things in the "de_equation_of_state" section so
    # we supply them there.
    block[names.de_equation_of_state, "a"] = a
    block[names.de_equation_of_state, "w"] = w
    block[names.de_equation_of_state, "z"] = z


    # This indicates success; non-zero values mean that these parameters should have
    # zero likelihood and be rejected.
    return 0