from cosmosis.datablock import option_section, names
from scipy.interpolate import interp2d
import numpy as np
import baccoemu

def setup(options):
    mode = options.get_string(option_section, "mode", default='pk_nl_only')
    
    # do this only once in the setup phase
    emulator = baccoemu.Matter_powerspectrum()

    allowed_modes = ['pk_nl_only','pk_nl_plus_baryons']
    if mode not in allowed_modes:
        raise ValueError('unrecognised setting: %s'%mode)
    
    return mode,emulator

def check_params(params):
    if not (params['omega_cold']>0.15) & (params['omega_cold']<0.6):
        raise ValueError('Omega_c must be in range [0.15,0.6]')

    if not (params['omega_baryon']>0.03) & (params['omega_baryon']<0.07):
        raise ValueError('Omega_b must be in range [0.03,0.07]')

    if not (params['ns']>0.92) & (params['ns']<1.01):
        raise ValueError('ns must be in range [0.92,1.01]')

    if not (params['hubble']>0.6) & (params['hubble']<0.8):
        raise ValueError('h must be in range [0.6,0.8]')

    if not (params['neutrino_mass']>0.) & (params['neutrino_mass']<0.4):
        raise ValueError('mnu must be in range [0,0.4]')

    return 0

def execute(block, config):
    mode, emulator = config

    if mode=='pk_nl_only':
        run_emulator_no_baryons(block,emulator)
    elif mode=='pk_nl_plus_baryons':
        run_emulator_with_baryons(block,emulator)
    else:
        raise ValueError('This should not happen')

    return 0

def run_emulator_with_baryons(block,emulator):

    P_lin = block['matter_power_lin','p_k']
    k_lin = block['matter_power_lin','k_h']
    z_lin = block['matter_power_lin','z']
  
    # This is required to avoid a bounds error 
    zmask = z_lin<1.5
    a_lin = 1/(1+z_lin[zmask])

    kmask = (k_lin<4.9) & (k_lin>0.0001)

    params = {
        'omega_cold'    :  block['cosmological_parameters','omega_m'] - block['cosmological_parameters','omega_nu'],
        'A_s'           :  block['cosmological_parameters','a_s'],
        'omega_baryon'  :  block['cosmological_parameters','omega_b'],
        'ns'            :  block['cosmological_parameters','n_s'],
        'hubble'        :  block['cosmological_parameters','h0'],
        'neutrino_mass' :  block['cosmological_parameters','mnu'],
        'w0'            :  block['cosmological_parameters','w'],
        'wa'            :  0.0,
        'expfactor'     :  a_lin,
        'M_c'           :  block['baryon_parameters','m_c'],
        'eta'           :  block['baryon_parameters','eta'],
        'beta'          :  block['baryon_parameters','beta'],
        'M1_z0_cen'     :  block['baryon_parameters','m1_z0_cen'],
        'theta_out'     :  block['baryon_parameters','theta_out'],
        'theta_inn'     :  block['baryon_parameters','theta_inn'],
        'M_inn'         :  block['baryon_parameters','m_inn']}

    # check we're within the allowed bounds for the emulator
    check_params(params)
    
    k_nl = k_lin
    z_nl = z_lin

    # evaluate the nonlinear growth factor as a function of k and z
    k_bacco, F = emulator.get_nonlinear_boost(k=k_lin[kmask], cold=False, **params)
    
    I_nl = interp2d(np.log10(k_bacco), z_lin[zmask], np.log10(F))
    F_interp = 10**I_nl(np.log10(k_nl),z_nl)

    # same thing for baryons
    k_bacco, S = emulator.get_baryonic_boost(k=k_lin[kmask], **params)

    I_baryon = interp2d(np.log10(k_bacco), z_lin[zmask], S)
    S_interp = I_baryon(np.log10(k_nl),z_nl)
    
    # apply the factor
    P_nl = F_interp * S_interp * P_lin

    block._copy_section('matter_power_lin','matter_power_nl')

    # save everything
    block['matter_power_nl','p_k'] = P_nl
    block['matter_power_nl','k_h'] = k_nl
    block['matter_power_nl','k'] = z_nl

    return 0


    
def run_emulator_no_baryons(block,emulator):

    P_lin = block['matter_power_lin','p_k']
    k_lin = block['matter_power_lin','k_h']
    z_lin = block['matter_power_lin','z']

    # bacco specific stuff
    # This is required to avoid a bounds error
    zmask = z_lin<1.5
    a_lin = 1/(1+z_lin[zmask])

    kmask = (k_lin<4.9) & (k_lin>0.0001)

    Omm = block['cosmological_parameters','omega_m']
    Omb = block['cosmological_parameters','omega_b']
    Omnu = block['cosmological_parameters','omega_nu']

    params = {
        'omega_cold'    :  Omm-Omnu,
        'A_s'   :  block['cosmological_parameters','a_s'],
        'omega_baryon'  :  Omb,
        'ns'            :  block['cosmological_parameters','n_s'],
        'hubble'        :  block['cosmological_parameters','h0'],
        'neutrino_mass' :  block['cosmological_parameters','mnu'],
        'w0'            : block['cosmological_parameters','w'],
        'wa'            :  0.0,
        'expfactor'     :  a_lin}

    # check we're within the allowed bounds for the emulator
    check_params(params)

    # evaluate the nonlinear growth factor as a function of k and z
    k_bacco, F = emulator.get_nonlinear_boost(k=k_lin[kmask], cold=False, **params)

    k_nl = k_lin
    z_nl = z_lin

    # interpolate it to the same sampling as the linear matter power spectrum
    I = interp2d(np.log10(k_bacco), z_lin[zmask], np.log10(F))
    F_interp = 10**I(np.log10(k_nl),z_nl)

    # apply the factor
    P_nl = F_interp * P_lin


    block._copy_section('matter_power_lin','matter_power_nl')

    # save everything
    block['matter_power_nl','p_k'] = P_nl
    block['matter_power_nl','k_h'] = k_nl
    block['matter_power_nl','k'] = z_nl

    return 0
