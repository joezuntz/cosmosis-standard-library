"""
An interface to the baryonic code.

"""
from cosmosis.datablock import option_section, names as section_names
import pyspk.model as spk
import numpy as np
import astropy.cosmology
import warnings
import sys
from scipy.interpolate import LinearNDInterpolator


def my_ceil(a, precision=0):
    return np.true_divide(np.ceil(a * 10**precision), 10**precision)


def my_floor(a, precision=0):
    return np.true_divide(np.floor(a * 10**precision), 10**precision)


round_up = np.vectorize(my_ceil)
round_down = np.vectorize(my_floor)


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

    try:
        fb_table = options[option_section, 'fb_table']
    except:
        fb_table = None

    extrapolate = options.get_bool(option_section, 'extrapolate', default=False)

    # These are parameters for spk
    config = {}
    config["verbose"] = verbose
    config["SO"] = SO
    config["fb_table"] = fb_table
    config["extrapolate"] = extrapolate


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
            raise ValueError(f"[SPK]: Unknown astropy model: {model}")        

    config["model_class"] = model_class

    return config

def execute(block, config):
    verbose = config['verbose']
    SO = config['SO']
    fb_table = config['fb_table']
    extrapolate = config['extrapolate']

    # Create our cosmological model
    model_class = config["model_class"]

    M_halo = None
    fb = None
    cosmo = None
    if model_class:
            cosmo = set_params(block, model_class)
        
    if not verbose:
        warnings.filterwarnings("ignore", category=UserWarning)

    # Load the current non-linear matter power spectrum
    section = section_names.matter_power_nl
    k, z_array, P = (block.get_grid(section, "k_h", "z", "P_k"))

    params = ['fb_a', 'fb_pow', 'fb_pivot', 'alpha', 'beta', 'gamma']

    spk_params = {}
    for param_i in params:
        try:
            spk_params[param_i] = block.get('spk', param_i)
        except:
            spk_params[param_i] = None

    modes = [
        ('Power law', ['fb_a', 'fb_pow', 'fb_pivot'], spk_params['fb_a'], spk_params['fb_pow'], spk_params['fb_pivot']),
        ('Non-parametric fb table', ['fb_table'], fb_table),
        ('Redshift dependent power law', ['alpha', 'beta', 'gamma', 'astropy_model'], spk_params['alpha'], spk_params['beta'], spk_params['gamma'], cosmo)
    ]
    provided_modes = sum(any(param is not None for param in mode[2:]) for mode in modes)

    if provided_modes != 1:
        provided_params = [param for mode in modes for param in mode[1] if mode[2:].count(None) != len(mode[2:])]
        raise ValueError(f"[SPK]: Only one mode to specify the baryon fraction should be provided. You provided: {provided_params}")

    for mode in modes:
        if any(param is not None for param in mode[2:]):
            missing_params = [param_name for param_name, param_value in zip(mode[1], mode[2:]) if param_value is None]
            if len(missing_params) != 0:
                raise ValueError(f"[SPK]: The following parameter(s) is(are) missing in mode '{mode[0]}': {', '.join(missing_params)}.")

    if fb_table:
        tab = np.loadtxt(fb_table, skiprows=1, delimiter=',')
        inter = LinearNDInterpolator(tab[:, [0, 1]], tab[:, 2], rescale=True)

    sup_array = np.empty_like(P)

    for i in range(len(z_array)):
        z = z_array[i]

        if fb_table:
            min_mass = round_down(np.log10(spk.optimal_mass(SO, z, np.max(k))), 1)
            max_mass = round_up(np.log10(spk.optimal_mass(SO, z,  np.min(k))), 1)
            M_halo = np.logspace(min_mass, max_mass, 100)
            fb = inter(z, M_halo) 
            if np.isnan(np.sum(fb)):
                raise ValueError(f"[SPK]: Requested values (z and/or M_halo) outside of the convex hull of the input points in fb_table. Hint: Check your fb_table and the requested redshifts.")

        k, sup = spk.sup_model(SO=SO, z=z, fb_a=spk_params['fb_a'], fb_pow=spk_params['fb_pow'], 
                               fb_pivot=spk_params['fb_pivot'], M_halo=M_halo, 
                               fb=fb, extrapolate=extrapolate, 
                               alpha=spk_params['alpha'], beta=spk_params['beta'], 
                               gamma=spk_params['gamma'], cosmo=cosmo, k_array=k, verbose=verbose)
        sup_array[:, i] = sup

    P_mod = P * sup_array

    block.replace_double_array_nd(section, "P_k", P_mod.T)

    # All is good - return
    return 0
