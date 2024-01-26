import numpy as np
import copy
import pickle
from .utils import _transform_space, MyProgressBar

__all__ = ["get_baryon_fractions"]

def load_baryonic_emu(fold_name=None, detail_name=None, verbose=True):
    """Loads in memory the baryonic boost emulator, described in Aricò et al. 2020c.

    :param verbose: whether to activate the verbose mode, defaults to True.
    :type verbose: boolean, optional

    :return: a dictionary containing the emulator object
    :rtype: dict
    """
    import os
    if verbose:
        print('Loading Baryonic Emulator...')
    basefold = os.path.dirname(os.path.abspath(__file__))

    from tensorflow.keras.models import load_model

    old_emulator_names = [(basefold + '/' +
                         "NN_emulator_sg_0.99_15000_PCA6_BNFalse_DO0.0_NL2_data_midmfcorr_10k")]
    for old_emulator_name in old_emulator_names:
        if os.path.exists(old_emulator_name):
            import shutil
            shutil.rmtree(old_emulator_name)

    if fold_name is not None:
        emulator_name = fold_name
    else:
        emulator_name = (basefold + '/' +
                         "baryon_emu_1.0.0")

        if (not os.path.exists(emulator_name)):
            import urllib.request
            import tarfile
            import ssl
            ssl._create_default_https_context = ssl._create_unverified_context
            print('Downloading Emulator data (5 Mb)...')
            urllib.request.urlretrieve(
                'https://bacco.dipc.org/baryon_emu_1.0.0.tar',
                emulator_name + '.tar',
                MyProgressBar())
            tf = tarfile.open(emulator_name+'.tar', 'r')
            tf.extractall(path=basefold)
            tf.close()
            os.remove(emulator_name + '.tar')

    emulator = {}
    emulator['emu_type'] = 'nn'
    emulator['model'] = load_model(emulator_name, compile=False)
    detail_name = 'k_scaler_bounds.pkl' if detail_name is None else detail_name
    with open(emulator_name + '/' +detail_name, 'rb') as f:
        emulator['scaler'] = pickle.load(f)
        emulator['pca'] = pickle.load(f)
        emulator['k'] = pickle.load(f)
        emulator['values'] = pickle.load(f)
        if fold_name == None:
            emulator['rotation'] = pickle.load(f)
        emulator['bounds'] = pickle.load(f)

    if fold_name == None:
        emulator['keys'] = ['omega_cold', 'sigma8_cold', 'omega_baryon',
             'ns', 'hubble', 'neutrino_mass', 'w0', 'wa', 'M_c', 'eta', 'beta', 'M1_z0_cen',
              'theta_inn', 'M_inn', 'theta_out', 'expfactor']
    else:
        emulator['keys'] = ['omega_cold', 'sigma8_cold', 'omega_baryon',
             'hubble', 'M_c', 'eta', 'beta', 'M1_z0_cen',
              'theta_inn','M_inn','theta_out', 'expfactor']
    if verbose:
        print('Baryonic Emulator loaded in memory.')
    return emulator

def get_baryon_fractions(M_200c, omega_cold=None, omega_matter=None, omega_baryon=None,
                         sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
                         w0=None, wa=None, expfactor=None, M_c=None, eta=None, beta=None,
                         M1_z0_cen=None, theta_out=None, theta_inn=None, M_inn=None):
    """Compute the mass fraction of the different baryonic components, following the
    baryonic correction model described in Aricò et al 2020b (see also Aricò et al 2020a).

    :param M_200: Halo mass inside a sphere which encompass a density 200 times larger than the critical density of the Universe.
    :type array_like
    :param omega_cold: omega cold matter (cdm + baryons), either omega_cold
                        or omega_matter should be specified, if both are specified
                        they should be consistent
    :type omega_cold: float or array
    :param omega_matter: omega total matter (cdm + baryons + neutrinos), either omega_cold
                        or omega_matter should be specified, if both are specified
                        they should be consistent
    :type omega_matter: float or array
    :param sigma8_cold: rms of cold (cdm + baryons) linear perturbations, either sigma8_cold
                        or A_s should be specified, if both are specified they should be
                        consistent
    :type sigma8_cold: float or array
    :param A_s: primordial scalar amplitude at k=0.05 1/Mpc, either sigma8_cold
                or A_s should be specified, if both are specified they should be
                consistent
    :type A_s: float or array
    :param hubble: adimensional Hubble parameters, h=H0/(100 km/s/Mpc)
    :type hubble: float or array
    :param ns: scalar spectral index
    :type ns: float or array
    :param neutrino_mass: total neutrino mass
    :type neutrino_mass: float or array
    :param w0: dark energy equation of state redshift 0 parameter
    :type w0: float or array
    :param wa: dark energy equation of state redshift dependent parameter
    :type wa: float or array
    :param expfactor: expansion factor a = 1 / (1 + z)
    :type expfactor: float or array
    :param M_c: mass fraction of hot gas in haloes
    :type M_c: float or array
    :param eta: extent of ejected gas
    :type eta: float or array
    :param beta: mass fraction of hot gas in haloes
    :type beta: float or array
    :param M1_z0_cen: characteristic halo mass scale for central galaxy
    :type M1_z0_cen: float or array
    :param theta_out: density profile of hot gas in haloes
    :type theta_out: float or array
    :param theta_inn: density profile of hot gas in haloes
    :type theta_inn: float or array
    :param M_inn: density profile of hot gas in haloes
    :type M_inn: float or array

    :return: a dictionary containing the baryonic components mass fractions, with the following keywords:

    #. 'ej_gas' -> ejected gas
    #. 'cen_galaxy' -> central galaxy
    #. 'sat_galaxy' -> satellite galaxy
    #. 'bo_gas' -> bound gas
    #. 're_gas' -> reaccreted gas
    #. 'dark_matter' -> dark matter
    #. 'gas' -> total gas
    #. 'stellar' -> total stars
    #. 'baryon' -> total baryons

    :rtype: dict
    """
    _kwargs = locals()
    kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}

    coord = copy.deepcopy(kwargs)
    for par in ['M_c','eta','beta','M1_z0_cen','theta_out','theta_inn','M_inn']:
        coord[par] = 10**(coord[par])

    stellar_pars = _SMHM_relation(coord)
    frac_cg = _galaxy_fraction(stellar_pars, M_200c, 'centrals')
    frac_sg = _galaxy_fraction(stellar_pars, M_200c, 'satellites')
    frac_stars = frac_cg + frac_sg


    o_m = coord['omega_matter'] if coord['omega_matter'] is not None else coord['omega_cold'] + coord['neutrino_mass'] / 93.14 / coord['hubble']**2

    baryon_matter_ratio = coord['omega_baryon']/o_m

    if np.any(frac_stars > baryon_matter_ratio):
        raise ValueError(""" Your stellar fraction is larger than the baryon fraction!
                              Please use meaningful stellar parameters""")

    frac_bg = (baryon_matter_ratio-frac_stars)/(1+(coord['M_c']/M_200c)**coord['beta'])
    frac_bar = frac_bg + frac_stars

    M_r = 1e16
    beta_r = 2.
    frac_re = (baryon_matter_ratio-frac_bar)/(1+(M_r/M_200c)**beta_r)
    frac_bg -= frac_re #fraction of reaccreted gas is taken from the bound gas
    frac_eg = baryon_matter_ratio - frac_bar
    frac_dm = 1 - baryon_matter_ratio

    frac_gas = frac_bg + frac_eg + frac_re
    return {'ej_gas':      np.array(frac_eg, dtype=np.float32),
            'cen_galaxy':  np.array(frac_cg, dtype=np.float32),
            'sat_galaxy':  np.array(frac_sg, dtype=np.float32),
            'bo_gas':      np.array(frac_bg, dtype=np.float32),
            're_gas':      np.array(frac_re, dtype=np.float32),
            'dark_matter': np.array(frac_dm, dtype=np.float32),
            'gas':  np.array(frac_gas, dtype=np.float32),
            'stellar':  np.array(frac_stars, dtype=np.float32),
            'baryon':  np.array(baryon_matter_ratio, dtype=np.float32),
            }

def _SMHM_relation(coordinates):
    """
    Internal function which evolve the Stellar Mass to Halo Mass (SMHM) relation parameters
    at the correct redshift, following Behroozi et al. 2013.

    :param coordinates: a set of coordinates in parameter space
    :type coordinates: dict
    :param a: expansion factor
    :type a: float
    :return: dictionary with the evolved SHAM parameters.
    :rtype: dictionary
    """
    a = coordinates['expfactor']
    z = 1/a -1
    nu = np.exp(-4*a**2) #Exponential cutoff of evolution of M ∗ (M h ) with scale factor

    pars = {}
    pars['M1_z0_cen'] = np.float64(coordinates['M1_z0_cen'])
    pars['epsilon_z0_cen'] = np.float32(0.023)
    pars['alpha_z0_cen'] = np.float32(-1.779)
    pars['gamma_z0_cen'] = np.float32(0.547)
    pars['delta_z0_cen'] = np.float32(4.394)

    pars['M1_fsat'] =       np.float32(1.59)
    pars['epsilon_fsat'] =  np.float32(0.2)
    pars['alpha_fsat'] =    np.float32(0.16)
    pars['gamma_fsat'] =    np.float32(1.67)
    pars['delta_fsat'] =    np.float32(0.99)

    ini = {'a':0,'a2':0,'z':0}
    stellar_pars = {
                'M1': copy.deepcopy(ini),      # Charactheristic halo mass Msolar/h
                'epsilon':copy.deepcopy(ini),   # Characteristic stellar mass to halo mass ratio
                'alpha':copy.deepcopy(ini),    # Faint-end slope of SMHM relation
                'delta':copy.deepcopy(ini),    #Index of subpower law at massive end of SMHM relation
                'gamma':copy.deepcopy(ini)    #Strength of subpower law at massive end of SMHM relation
                }

    stellar_pars['M1']['a'] = -1.793;  stellar_pars['M1']['z'] = -0.251
    stellar_pars['alpha']['a'] = 0.731;
    stellar_pars['gamma']['a'] = 1.319; stellar_pars['gamma']['z'] = 0.279
    stellar_pars['delta']['a'] = 2.608; stellar_pars['delta']['a'] = -0.043
    stellar_pars['epsilon']['a'] = -0.006;  stellar_pars['epsilon']['a2'] = -0.119;

    #for central galaxies:
    for p in stellar_pars.keys():
        p0_cen = np.log10(pars[p+'_z0_cen']) if p in ['M1','epsilon'] else pars[p+'_z0_cen']
        pcen_z =  p0_cen + nu*(stellar_pars[p]['a']*(a-1) + stellar_pars[p]['z']*z) + stellar_pars[p]['a2']*(a-1)
        pars[p+'_cen'] = 10**(pcen_z) if p in ['M1','epsilon'] else pcen_z

    #for satellites:
    for p in stellar_pars.keys():
        pars[p+'_sat'] = pars[p+'_cen'] * pars[p+'_fsat']
    return pars

def _galaxy_fraction(pars, M_200, type='centrals'):
    """
    Function which compute the galaxy mass fractions.
    :param pars: set of SHAM parameters
    :type pars: dict
    :param M_200: Halo mass inside a sphere which encompass a density 200 times larger than the critical density of the Universe.
    :type array_like
    :param type: 'centrals' to select central galaxies, 'satellites' to select satellites galaxies.
    :type string
    :return: galaxy mass fractions.
    :rtype: array_like
    """
    t = '_cen' if type == 'centrals' else '_sat'
    return _stellar_fraction(M_200, M1=pars['M1'+t], alpha=pars['alpha'+t],
            gamma=pars['gamma'+t],delta=pars['delta'+t],epsilon=pars['epsilon'+t])

def _stellar_fraction(M_200, M1=1.526e11, alpha= -1.779, gamma = 0.547, delta = 4.394, epsilon=0.023):
        """
        Internal function which compute the galaxy mass fractions.

        :param M_200: Halo mass inside a sphere which encompass a density 200 times larger than the critical density of the Universe.
        :type array_like
        :return: galaxy mass fractions.
        :rtype: array_like
        """
        fact = 10**( _g(np.log10(M_200/M1), alpha, gamma, delta) - _g(0, alpha, gamma, delta))
        return epsilon * (M1/M_200) * fact

def _g(x, alpha= -1.779, gamma = 0.547, delta = 4.394):
    return -np.log10( 10**(alpha*x) + 1) + \
        delta * ( (np.log10( 1+np.exp(x) ))**gamma ) / ( 1 + np.exp(10**(-x)) )
