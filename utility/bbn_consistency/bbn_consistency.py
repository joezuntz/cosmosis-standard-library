from cosmosis.datablock import names, option_section
import numpy as np
from scipy import interpolate
import os

cosmo = names.cosmological_parameters
dirname = os.path.split(__file__)[0]


def setup(options):
    default_datafile = os.path.join(dirname, "helium.dat")
    # File name for BBN data
    datafile = options.get_string(
        option_section, "data", default=default_datafile)
    # Load and spline data
    dat = np.genfromtxt(datafile, names=True, comments='#')
    spline = interpolate.bisplrep(dat['ombh2'], dat['DeltaN'], dat['Yp'])

    input_name = options.get_string(
        option_section, "input_name", default="delta_neff")
    assert (input_name == "delta_neff") or (input_name == "massless_nu")

    # Saved data will be returned later
    return {'spline': spline, 'ombh2_min': np.min(dat['ombh2']), 'ombh2_max': np.max(dat['ombh2']), 'DeltaN_min': np.min(dat['DeltaN']), 'DeltaN_max': np.max(dat['DeltaN']), 
        'input_name':input_name}


def execute(block, t):
    input_name = t['input_name']
    # This loads a value from the section "cosmological_parameters" that we read above.
    ombh2 = block[cosmo, "ombh2"]
    if input_name == "delta_neff":
        if block.has_value(cosmo, "delta_neff"):
            DeltaN = block[cosmo, "delta_neff"]
        else:
            DeltaN = 0.0
    elif input_name == "massless_nu":
        DeltaN = block[cosmo, "massless_nu"] + block[cosmo, "massive_nu"] - 3.046


    # Check for out-of-range parameters
    if ombh2 < t['ombh2_min'] or ombh2 > t['ombh2_max'] or DeltaN < t['DeltaN_min'] or DeltaN > t['DeltaN_max']:
        return 1

    # save the helium fraction back to the block
    block[cosmo, "yhe"] = interpolate.bisplev(ombh2, DeltaN, t['spline'])

    # We tell CosmoSIS that everything went fine by returning zero
    return 0


def cleanup(config):
    # Python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
