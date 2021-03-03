from cosmosis.datablock import names, option_section
cosmo = names.cosmological_parameters

"""
Skip point in parameter space if w0+wa>=0  to avoid CAMB errors. 
"""

def setup(options):
    verbose = options.get_bool(option_section, "verbose", default=False)
    return {'verbose':verbose}

def execute(block, config):

    w0 = block[cosmo, "w"]
    wa = block[cosmo, "wa"]
    
    if w0 + wa >= 0:
        if config['verbose']:
            print("Parameters fail w0+wa<0 requirement: w0 = {}  wa = {}".format(w0, wa))
        return 1

    return 0


def cleanup(config):
    pass
