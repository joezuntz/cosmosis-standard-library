from cosmosis.datablock import option_section, names
import numpy as np
import os

dirname = os.path.split(__file__)[0]
default_data_dir = os.path.join(dirname, "data")

def setup(options):
    data_dir = options.get_string(option_section, "data_dir", default_data_dir)
    replace_pp = options.get_bool(option_section, "replace_pp", False)
    zero_remainder = options.get_bool(option_section, "zero_remainder", False)

    spectra = ["TT", "TE", "EE", "BB"]
    
    if replace_pp:
        spectra.append("PP")

    config = {
        'spectra' : spectra,
        'zero_remainder' : zero_remainder,
    }

    for name in spectra:
        filename = os.path.join(data_dir, 'cl{}'.format(name))
        cl = np.loadtxt(filename)
        ell = np.arange(len(cl))

        # Cosmosis stores these as l(l+1)cl/2pi
        cl *= ell*(ell+1.)/2./np.pi

        # Cosmosis does not save the monopole or dipole
        cl = cl[2:]
        print("Loaded fiducial {} from {}".format(name, filename))
        print("Will replace first {} values of TT spectra".format(len(cl)))
        config[name] = cl

    return config


def execute(block, config):
    spectra = config['spectra']
    for name in spectra:

        # Replace the first chunk of the spectrum
        # with the version loaded from file.
        cl = config[name]
        current_cl = block[names.cmb_cl, name]
        n_replaced = cl.size
        current_cl[:n_replaced] = cl

        # Optionally zero out the rest of the spectrum, where our
        # fiducial values do not go to high enough ell.
        # Otherwise use the version calculated in camb.
        if config['zero_remainder']:
            current_cl[n_replaced:] = 0.0

        # Store the result
        block[names.cmb_cl, name] = current_cl
    return 0