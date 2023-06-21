"""
An interface to the baryonic code.

"""
from cosmosis.datablock import option_section, names as section_names
import pyspk.model as spk
import numpy as np
import astropy.cosmology
import warnings
import sys


def set_params(block, model_class):
    # These pairs are the astropy parameter name and the cosmosis param name
    params_needed_by_class = {  
        astropy.cosmology.FlatLambdaCDM: [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b')],


        astropy.cosmology.FlatwCDM:      [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), 
                                          ('w0', 'w')],

        astropy.cosmology.Flatw0waCDM:   [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), 
                                          ('w0', 'w'), ('wa', 'wa')],

        astropy.cosmology.LambdaCDM:     [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), ('Ode0', 'omega_lambda')],


        astropy.cosmology.wCDM:          [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), ('Ode0', 'omega_lambda'),
                                          ('w0', 'w')],

        astropy.cosmology.w0waCDM:       [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), ('Ode0', 'omega_lambda'), 
                                          ('w0', 'w'), ('wa', 'wa')],

        astropy.cosmology.w0wzCDM:       [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), ('Ode0', 'omega_lambda'),
                                          ('w0', 'w'), ('wz', 'wz')],

        astropy.cosmology.w0wzCDM:       [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), ('Ode0', 'omega_lambda'),
                                          ('w0', 'w'), ('wz', 'wz')],

        astropy.cosmology.wpwaCDM:       [('H0', 'hubble'), ('Om0', 'omega_m'), ('Ob0', 'omega_b'), ('Ode0', 'omega_lambda'), 
                                          ('wp', 'wp'),('wa', 'wa'), ('zp', 'zp')],

    }

    params_needed = params_needed_by_class[model_class]

    # Pull the parameters we need out from cosmosis
    params = {}
    for astropy_name, cosmosis_name in params_needed:
        params[astropy_name] = block[section_names.cosmological_parameters, cosmosis_name]
    
    # Create the astropy object that does the calculations
    cosmo = model_class(**params)

    return cosmo


def setup(options):
    verbose = options.get_bool(option_section, 'verbose', default=False)
    SO = options.get_int(option_section, 'SO', default=500)
    model = options.get_string(option_section, "astropy_model", default='None')

    # These are parameters for spk
    config = {}

    config["verbose"] = verbose
    config["SO"] = SO


    # Models available in astropy.cosmology
    astropy_models = {
        "flatlambdacdm": astropy.cosmology.FlatLambdaCDM,
        "flatw0wacdm": astropy.cosmology.Flatw0waCDM,
        "flatwcdm": astropy.cosmology.FlatwCDM,
        "lambdacdm": astropy.cosmology.LambdaCDM,
        "w0wacdm": astropy.cosmology.w0waCDM,
        "w0wzcdm": astropy.cosmology.w0wzCDM,
        "wcdm": astropy.cosmology.wCDM,
        "wpwacdm": astropy.cosmology.wpwaCDM,
    }

    # Find the model the user has specified
    model_class = None

    if 'None' not in model:
        model_class = astropy_models.get(model.lower())
        if model_class is None:
            raise ValueError("Unknown astropy model {}".format(model))        

    config["model_class"] = model_class

    return config

def execute(block, config):
    verbose = config['verbose']
    SO = config['SO']

    #Create our cosmological model
    
    model_class = config["model_class"]

    cosmo = None
    if model_class:
            cosmo = set_params(block, model_class)
        
    if not verbose:
        warnings.filterwarnings("ignore", category=UserWarning)

    # Load the current unmodulated matter power spectrum
    section = section_names.matter_power_nl
    k, z_array, P = (block.get_grid(section, "k_h", "z", "P_k"))

    params = ['fb_a', 'fb_pow', 'fb_pivot', 'extrapolate', 'alpha', 'beta', 'gamma', 'M_halo', 'fb']

    spk_params = {}
    for param in params:
        try:
            spk_params[param] = block.get('spk', param)
        except:
            spk_params[param] = None


    sup_array = np.empty_like(P)
    # print(np.shape(sup_array))
    for i in range(len(z_array)):
        z = z_array[i]

        k, sup = spk.sup_model(SO=SO, z=z, fb_a=spk_params['fb_a'], fb_pow=spk_params['fb_pow'], 
                               fb_pivot=spk_params['fb_pivot'], M_halo=spk_params['M_halo'], 
                               fb=spk_params['fb'], extrapolate=spk_params['extrapolate'], 
                               alpha=spk_params['alpha'], beta=spk_params['beta'], 
                               gamma=spk_params['gamma'], cosmo=cosmo, k_array=k, verbose=verbose)
        sup_array[:, i] = sup

    P_mod = P * sup_array

    block.replace_double_array_nd(section, "P_k", P_mod.T)

    # All is good - return
    return 0
