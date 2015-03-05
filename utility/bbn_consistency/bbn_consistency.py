from cosmosis.datablock import names, option_section
import numpy as np
from scipy import interpolate

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.

    #Use this syntax to get a single parameter from the ini file section
    #for this module.  There is no type checking here - you get whatever the user
    #put in.
    
    datafile = options[option_section, "data"]

    #The call above will crash if "mode" is not found in the ini file.
    #Sometimes you want a default if nothing is found:
    ##high_accuracy = options.get(option_section, "high_accuracy", default=False)

    #Now you have the input options you can do any useful preparation
    #you want.  Maybe load some data, or do a one-off calculation.
    dat = np.genfromtxt(datafile, names=True, comments='#')
    spline = interpolate.bisplrep(dat['ombh2'], dat['DeltaN'], dat['Yp'])

    #Whatever you return here will be saved by the system and the function below
    #will get it back.  You could return 0 if you won't need anything.
    return {'spline':spline, 'ombh2_min':np.min(dat['ombh2']), 'ombh2_max':np.max(dat['ombh2']), 'DeltaN_min':np.min(dat['DeltaN']), 'DeltaN_max':np.max(dat['DeltaN'])}


def execute(block, t):
    #This function is called every time you have a new sample of cosmological and other parameters.
    #It is the main workhorse of the code. The block contains the parameters and results of any 
    #earlier modules, and the config is what we loaded earlier.

    #This loads a value from the section "cosmological_parameters" that we read above.
    ombh2 = block[cosmo, "ombh2"]
    if block.has_value(cosmo, "delta_neff"):
        DeltaN = block[cosmo, "delta_neff"]
    else:
        DeltaN = 0.0

    if ombh2 < t['ombh2_min'] or ombh2 > t['ombh2_max'] or DeltaN < t['DeltaN_min'] or DeltaN > t['DeltaN_max']:
        return 1

    # save the result back to the block like this.
    block[cosmo, "yhe"] = interpolate.bisplev(ombh2, DeltaN, t['spline'])

    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
