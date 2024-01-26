import numpy as np
import copy
import pickle
import os
from .utils import _transform_space, MyProgressBar, mean_absolute_exp_percentage_error, accuracy_exp_01, accuracy_exp_005
from scipy import interpolate
import tensorflow
tensorflow.compat.v1.logging.set_verbosity(tensorflow.compat.v1.logging.ERROR)
from tensorflow.keras.models import load_model
gpus = tensorflow.config.experimental.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tensorflow.config.experimental.set_memory_growth(gpu, True)
from .matter_powerspectrum import load_smeared_bao_emu


__all__ = ["Lbias_expansion"]

class Lbias_expansion(object):
    """
    A class to load and call the baccoemu for the Lagrangian bias expansion terms.
    By default, the nonlinear Lagrangian bias expansion terms emulator (described
    in Zennaro et al, 2021) and LPT lagrangian bias expansion terms emulator
    (described in Aricò et al, 2021) are loaded.
    Another emulator loaded by deafult is the IR-resummed linear power spectrum.

    :param lpt: whether to load the LPT emulator, defaults to True
    :type lpt: boolean, optional
    :param compute_sigma8: whether to load the sigma8 emulator, defaults to True
    :type compute_sigma8: boolean, optional
    :param smeared_bao: whether to load the smeared bao, defaults to True
    :type smeared_bao: boolean, optional
    :param nonlinear_boost: whether to load the nonlinear boost emulator, defaults to True
    :type nonlinear_boost: boolean, optional
    :param compute_sigma8: whether to load the sigma8 emulator, defaults to True
    :type compute_sigma8: boolean, optional
    :param verbose: whether to activate the verbose mode, defaults to True
    :type verbose: boolean, optional

    """
    def __init__(self, lpt=True, smeared_bao=True, nonlinear_boost=True, compute_sigma8=True, verbose=True):

        self.verbose = verbose

        self.compute_lpt = True if lpt else False
        self.compute_smeared_bao = True if smeared_bao else False

        self.cosmo_keys = np.array(['omega_cold', 'sigma8_cold', 'omega_baryon', 'ns',
                                'hubble', 'neutrino_mass', 'w0', 'wa', 'expfactor'])

        self.lb_term_labels = [r'$1 1$', r'$1 \delta$', r'$1 \delta^2$', r'$1 s^2$',
                               r'$ 1 \nabla^2\delta$', r'$\delta \delta$',
                               r'$\delta \delta^2$', r'$\delta s^2$',
                               r'$\delta \nabla^2\delta$', r'$\delta^2 \delta^2$',
                               r'$\delta^2 s^2$', r'$\delta^2 \nabla^2\delta$',
                               r'$s^2 s^2$', r'$s^2 \nabla^2\delta$',
                               r'$\nabla^2\delta \nabla^2\delta$']

        self.emulator = {}
        if self.compute_lpt:
            self.emulator['lpt'] = load_lpt_emu()

        if self.compute_smeared_bao:
            self.emulator['smeared_bao'] = load_smeared_bao_emu()

        self.compute_nonlinear_boost = True if nonlinear_boost else False
        if self.compute_nonlinear_boost:
            self.emulator['nonlinear'] = load_nonlinear_lbias_emu()

        self.compute_sigma8 = True if compute_sigma8 else False

        if self.compute_sigma8:
            from .matter_powerspectrum import Matter_powerspectrum
            self.matter_powerspectrum_emulator = Matter_powerspectrum(linear=False, smeared_bao=False,
            nonlinear_boost=False, baryonic_boost=False,
            compute_sigma8=True, verbose=verbose)

    def _get_parameters(self, coordinates, which_emu, grid=None):
        """
        Function that returns a dictionary of cosmological parameters,
        computing derived cosmological parameters, if not
        already present in the given coordinates, and checking the relevant boundaries.
        :param coordinates: a set of coordinates in parameter space
        :type coordinates: dict
        :param which_emu: kind of emulator: options are 'linear', 'nonlinear','baryon','smeared_bao','sigma8'
        :type which_emu: str
        :param grid: dictionary with parameter and vector of values where to evaluate the emulator, defaults to None
        :type grid: array_like, optional
        :return: coordinates with derived parameters
        :rtype: dict
        """
        coordinates = {key: np.atleast_1d(coordinates[key]) for key in set(list(coordinates.keys())) - set(['k', 'k_lin', 'pk_lin'])}

        avail_pars = [coo for coo in coordinates.keys() if coordinates[coo][0] is not None] #parameters currently available
        eva_pars = self.emulator[which_emu]['keys']  #parameters strictly needed to evaluate the emulator
        req_pars = self.emulator[which_emu]['keys']  #parameters needed for a computation
        comp_pars = list(set(req_pars)-set(avail_pars)) #parameters to be computed
        deriv_pars = ['omega_cold','sigma8_cold', 'A_s'] #derived parameters that can be computed
        miss_pars = list(set(comp_pars)-set(deriv_pars)) #parameters missing from coordinates
        extra_pars = list(set(req_pars)-set(eva_pars)) #requested parameters not needed for evaluation
        if miss_pars:
            print(f"{which_emu} emulator:")
            print(f"  Please add the parameter(s) {miss_pars} to your coordinates!")
            raise KeyError(f"{which_emu} emulator: coordinates need the following parameters: ", miss_pars)

        if ('omega_cold' in avail_pars) & ('omega_matter' in avail_pars):
            assert len(coordinates['omega_cold']) == len(coordinates['omega_matter']), 'Both omega_cold and omega_matter were provided, but they have different len'
            om_from_oc = coordinates['omega_cold'] + coordinates['neutrino_mass'] / 93.14 /coordinates['hubble']**2
            assert np.all(np.abs(coordinates['omega_matter'] - om_from_oc) < 1e-4), 'Both omega_cold and omega_matter were provided, but they are inconsistent among each other'

        if 'omega_cold' in comp_pars:
            if 'omega_matter' not in avail_pars:
                raise KeyError('One parameter between omega_matter and omega_cold must be provided!')

            omega_nu = coordinates['neutrino_mass'] / 93.14 /coordinates['hubble']**2
            coordinates['omega_cold'] = coordinates['omega_matter'] - omega_nu

        if ('sigma8_cold' not in avail_pars) & ('A_s' not in avail_pars):
            raise KeyError('One parameter between sigma8_cold and A_s must be provided!')

        if ('sigma8_cold' in  avail_pars) & ('A_s' in avail_pars):
            #commented for the cases where one is computed and same value is repeated
            #assert len(np.atleast_1d(coordinates['sigma8_cold'])) == len(atleast_1d(coordinates['A_s'])), 'Both sigma8_cold and A_s were provided, but they have different len'

            ignore_s8_pars = copy.deepcopy(coordinates)
            del ignore_s8_pars['sigma8_cold']
            s8_from_A_s = self.matter_powerspectrum_emulator.get_sigma8(**ignore_s8_pars)
            assert np.all(np.abs(coordinates['sigma8_cold'] - s8_from_A_s) < 1e-4), 'Both sigma8_cold and A_s were provided, but they are inconsistent among each other'

        if 'sigma8_cold' in comp_pars:
            tmp_coords = copy.deepcopy(coordinates)
            tmp_coords['cold']=True
            coordinates['sigma8_cold'] = np.atleast_1d(self.matter_powerspectrum_emulator.get_sigma8(**tmp_coords))

        if 'A_s' in comp_pars:
            tmp_coords = copy.deepcopy(coordinates)
            del tmp_coords['sigma8_cold']
            tmp_coords['A_s'] = 2e-9
            tmp_coords['cold'] = True
            coordinates['A_s'] = np.atleast_1d((coordinates['sigma8_cold'] / self.matter_powerspectrum_emulator.get_sigma8(**tmp_coords))**2 * tmp_coords['A_s'])


        pp = np.squeeze([coordinates[p][0] for p in eva_pars])
        coords_out = copy.deepcopy(coordinates)

        grid = {}
        for key in coordinates.keys():
            if len(np.atleast_1d(coordinates[key])) > 1:
                grid[key] = np.array(coordinates[key])

        if len(list(grid.keys()))==0:
            grid = None
        else:
            grid_structure = []
            for key in grid.keys():
                grid_structure.append(len(grid[key]))
            grid_structure = np.array(grid_structure)
            values, counts = np.unique(grid_structure, return_counts=True)
            counts_but_highest = np.delete(counts, np.argmax(counts))
            assert np.all(counts == counts[0]) | np.all(counts_but_highest == 1), 'When passing multiple coordinate sets you should either vary only on parameter, or all parameters should have the same len'

        if grid is not None:
            grid_pars = list(grid.keys()) # list of parameters that are varyied in a grid
            N = len(grid[grid_pars[0]])
            pp = np.tile(pp, (N, 1))
            for par in grid_pars:
                if par in eva_pars:
                    index = eva_pars.index(par)
                    pp[:,index] = np.float64(grid[par])
                if par in req_pars:
                    coords_out[par] = grid[par]
            pp = np.float64(pp)

        for i,par in enumerate(eva_pars):
            val = pp[i] if grid is None else pp[:,i]
            message = 'Param {}={} out of bounds [{}, {}]'.format(
                par, val, self.emulator[which_emu]['bounds'][i][0], self.emulator[which_emu]['bounds'][i][1])

            assert np.all(val >= self.emulator[which_emu]['bounds'][i][0]) & np.all(val <= self.emulator[which_emu]['bounds'][i][1]), message

        if extra_pars:
            cc = np.squeeze([coords_out[p] for p in extra_pars])
            if None in cc:
                raise ValueError(f'None in parameters: {extra_pars} = {cc}!')

        return coords_out, pp, grid

    def get_galaxy_real_pk(self, bias=None, omega_cold=None, omega_matter=None, omega_baryon=None,
                          sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
                          w0=None, wa=None, expfactor=None, k=None, **kwargs):
        """Compute the predicted galaxy auto pk and galaxy-matter cross pk given a set of bias parameters

        :param bias: a list of bias parameters, including b1, b2, bs2, blaplacian
        :type bias: array-like
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
        :param k: a vector of wavemodes in h/Mpc at which the nonlinear boost will be computed, if None
                  the default wavemodes of the nonlinear emulator will be used, defaults to None
        :type k: array_like, optional
        :return: k and P(k), a list of the emulated 15 LPT Lagrangian bias expansion terms
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}

        import itertools

        assert len(bias) == 4, 'Please, pass a valid bias array, with b1, b2, bs2, blaplacian'

        k, pnn = self.get_nonlinear_pnn(**kwargs)
        bias = np.concatenate(([1], bias))
        prod = np.array(list(itertools.combinations_with_replacement(np.arange(len(bias)), r=2)))

        pgal_auto = 0
        for i in range(len(pnn)):
            fac = 2 if prod[i, 0] != prod[i,1] else 1
            pgal_auto += bias[prod[i, 0]] * bias[prod[i,1]] * fac * pnn[i]
        pgal_cross = np.dot(bias, pnn[:5])

        return k, pgal_auto, pgal_cross


    def get_nonlinear_pnn(self, omega_cold=None, omega_matter=None, omega_baryon=None,
                          sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
                          w0=None, wa=None, expfactor=None, k=None, **kwargs):
        """Compute the prediction of the nonlinear cold matter power spectrum.

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
        :param k: a vector of wavemodes in h/Mpc at which the nonlinear boost will be computed, if None
                  the default wavemodes of the nonlinear emulator will be used, defaults to None
        :type k: array_like, optional
        :return: k and P(k), a list of the emulated 15 LPT Lagrangian bias expansion terms
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}

        if not self.compute_nonlinear_boost:
            raise ValueError("Please enable the l-bias nonlinear boost!")

        coordinates, pp, grid = self._get_parameters(kwargs, 'nonlinear')
        emulator = self.emulator['nonlinear']

        n_log = [0,1,5,9,12,14]
        n_switch = [10 , 11]

        _pp = _transform_space(np.array([pp]), space_rotation=False, bounds=emulator['bounds'])

        tmp_coords = copy.deepcopy(kwargs)
        tmp_coords['k'] = emulator['k']
        _, pk_lpt = self.get_lpt_pk(**tmp_coords)
        _, pk_bao = self.get_smeared_bao_pk(**tmp_coords)

        pk_lpt[0] = pk_bao
        pk_lpt[1] = pk_bao
        pk_lpt[5] = pk_bao
        pk_lpt[4] = -pk_bao * emulator['k']**2
        pk_lpt[8] = -pk_bao * emulator['k']**2
        pk_lpt[14] = pk_bao * emulator['k']**4

        P_nn = []
        for n in range(15):
            prediction = emulator['model'][n](_pp.reshape(-1,9), training=False)
            # import pdb; pdb.set_trace()
            if n in n_switch:
                lpt_term = pk_lpt[n + 2]
            else:
                lpt_term = pk_lpt[n]
            if n in n_log:
                _this_P_nn = np.squeeze(np.exp(emulator['scaler'][n].inverse_transform(prediction)) * lpt_term)
            else:
                _this_P_nn = np.squeeze(emulator['scaler'][n].inverse_transform(prediction) * lpt_term)
            if n in n_switch:
                if len(pk_lpt.shape) == 2:
                    _this_P_nn[:25] = _this_P_nn[:25] * pk_lpt[n][:25] / lpt_term[:25]
                else:
                    _this_P_nn[:,:25] = _this_P_nn[:,:25] * pk_lpt[n,:,:25] / lpt_term[:,:25]
            P_nn.append(_this_P_nn)

        if k is not None:
            if max(k) > 0.75:
                raise ValueError(f"""
            The maximum k of the l-bias nonlinear emulator must be 0.75 h/Mpc:
            the current value is {max(k)} h/Mpc""")
            if (min(k) <= 1e-2)&(self.verbose):
                print("WARNING: the nonlinear emulator is extrapolating to k < 0.01 h/Mpc!")

            new_P_nn = []
            for n in range(15):
                unexpected_negative = np.any(P_nn[n] <= 0.0) # can happen when allowing extrapolation
                if (n in n_log) & (unexpected_negative is False):
                    new_P_nn.append(np.exp(interpolate.interp1d(np.log(emulator['k']), np.log(P_nn[n]), kind='cubic', axis=0 if grid is None else 1, fill_value='extrapolate')(np.log(k))))
                else:
                    new_P_nn.append(interpolate.interp1d(np.log(emulator['k']), P_nn[n], kind='cubic', axis=0 if grid is None else 1, fill_value='extrapolate')(np.log(k)))
            P_nn = np.array(new_P_nn)
        else :
            k = emulator['k']

        return k, P_nn

    def get_lpt_pk(self, omega_cold=None, omega_matter=None, omega_baryon=None,
                   sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
                   w0=None, wa=None, expfactor=None, k=None, **kwargs):
        """Compute the prediction of the 15 LPT Lagrangian bias expansion terms.


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
        :param k: a vector of wavemodes in h/Mpc at which the nonlinear boost will be computed, if None
                the default wavemodes of the nonlinear emulator will be used, defaults to None
        :type k: array_like, optional
        :return: k and P(k), a list of the emulated 15 LPT Lagrangian bias expansion terms
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}

        if not self.compute_lpt:
            raise ValueError("Please enable the lpt emulator!")

        emulator = self.emulator['lpt']
        coordinates, pp, grid = self._get_parameters(kwargs, 'lpt')

        sub = emulator['sub']
        scaler = emulator['scaler']

        P_nn = []
        for n in range(15):
            pred = emulator['model'][n](pp.reshape(-1,9), training=False)
            prediction = np.squeeze(scaler[n].inverse_transform(pred))
            P_nn.append(prediction)

        if k is not None:
            if max(k) > 0.75:
                raise ValueError(f"""
            The maximum k of the l-bias lpt emulator must be 0.75 h/Mpc:
            the current value is {max(k)} h/Mpc""")
            if (min(k) <= 1e-2)&(self.verbose):
                print("WARNING: the l-bias lpt emulator is extrapolating to k < 0.01 h/Mpc!")

            for n in range(15):
                p_interp = interpolate.interp1d(np.log(emulator['k']), P_nn[n], kind='cubic', axis=0 if grid is None else 1, fill_value='extrapolate',
                assume_sorted=True)
                P_nn[n] = p_interp(np.log(k))
        else :
            k = emulator['k']

        P_nn = np.array([np.exp(P_nn[n])-sub[n] for n in range(15)])
        return k, P_nn

    def get_smeared_bao_pk(self, omega_cold=None, omega_matter=None, omega_baryon=None,
                           sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
                           w0=None, wa=None, expfactor=None, k=None, **kwargs):
        """Evaluate the smeared bao emulator at a set of coordinates in parameter space.

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
        :param k: a vector of wavemodes in h/Mpc at which the nonlinear boost will be computed, if None
                the default wavemodes of the nonlinear emulator will be used, defaults to None
        :type k: array_like, optional
        :return: k and P(k), a list of the emulated 15 LPT Lagrangian bias expansion terms
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}

        if not self.compute_smeared_bao:
            raise ValueError("Please enable the smeared bao emulator!")

        emulator = self.emulator['smeared_bao']
        coordinates, pp, grid = self._get_parameters(kwargs, 'smeared_bao')

        ypred = emulator['model'](pp.reshape(-1,9), training=False)
        pk_bao = np.squeeze(np.exp(emulator['scaler'].inverse_transform(ypred)))
        if k is not None:
            if (max(k) > 30.)|(min(k) < 1e-3):
                raise ValueError(f"""
                    A minimum k > 0.001 h/Mpc and a maximum k < 30 h/Mpc
                    are required for the linear emulator:
                    the current values are {min(k)} h/Mpc and {max(k)} h/Mpc
                    """)

            else:
                pk_bao = np.exp(interpolate.interp1d(np.log(emulator['k']),np.log(pk_bao), kind='cubic', axis=0 if grid is None else 1, fill_value='extrapolate')(np.log(k)))
        else:
            k = emulator['k']
        return  k, pk_bao

    def _vector_eval_lpt_pk(self, coordinates):
        """Get the lpt pk in a fast vectorized way

        Note: no checks are performed to increase speed, be aware of what you do

        :param coordinates: 2D array with a grid of cooordinates
        :type coordinates: array
        """
        emulator = self.lpt_emulator
        sub = emulator['sub']
        scaler = emulator['scaler']

        P_nn = []
        for n in range(15):
            pred = emulator['model'][n](coordinates, training=False)
            prediction = np.exp(scaler[n].inverse_transform(pred)) - sub[n]
            P_nn.append(prediction)
        return np.array(P_nn)

    def _vector_eval_smeared_bao_pk(self, coordinates):
        """Get the nonlinear pk in a fast vectorized way

        Note: no checks are performed to increase speed, be aware of what you do

        :param coordinates: 2D array with a grid of cooordinates
        :type coordinates: array
        """
        emulator = self.smeared_bao_emulator
        ypred = emulator['model'](coordinates, training=False)
        pk_bao = np.exp(emulator['scaler'].inverse_transform(ypred))
        return  pk_bao


    def _vector_eval_nonlinear_pnn(self, coordinates):
        """Get the nonlinear pk in a fast vectorized way

        Note: no checks are performed to increase speed, be aware of what you do

        :param coordinates: 2D array with a grid of cooordinates
        :type coordinates: array
        """
        from scipy.interpolate import interp1d
        emulator = self.nonlinear_emulator
        _pp = _transform_space(coordinates, space_rotation=False, bounds=emulator['bounds'])

        pk_lpt = self._vector_eval_lpt_pk(coordinates)
        pk_bao = self._vector_eval_smeared_bao_pk(coordinates)
        pk_bao_interp = np.exp(interp1d(np.log(self.smeared_bao_emulator['k']), np.log(pk_bao), kind='cubic')(np.log(emulator['k'])))

        pk_lpt[0,:,:] = pk_bao_interp
        pk_lpt[1,:,:] = pk_bao_interp
        pk_lpt[5,:,:] = pk_bao_interp
        pk_lpt[4,:,:] = -pk_bao_interp * emulator['k']**2
        pk_lpt[8,:,:] = -pk_bao_interp * emulator['k']**2
        pk_lpt[14,:,:] = pk_bao_interp * emulator['k']**4

        P_nn = []
        n_log = [0,1,5,9,12,14]
        n_switch = [10 , 11]
        for n in range(15):
            prediction = emulator['model'][n](_pp, training=False)
            if n in n_switch:
                lpt_term = pk_lpt[n + 2,:,:]
            else:
                lpt_term = pk_lpt[n,:,:]
            if n in n_log:
                _this_P_nn = np.exp(emulator['scaler'][n].inverse_transform(prediction)) * lpt_term
            else:
                _this_P_nn = emulator['scaler'][n].inverse_transform(prediction) * lpt_term
            if n in n_switch:
                _this_P_nn[:,:25] = _this_P_nn[:,:25] * pk_lpt[n,:,:25] / lpt_term[:,:25]
            P_nn.append(_this_P_nn)

        return np.array(P_nn)

    def _vector_eval_galaxy_real_pk(self, coordinates, bias):
        """Get the galaxy pk in a fast vectorized way

        Note: no checks are performed to increase speed, be aware of what you do

        :param coordinates: 2D array with a grid of cooordinates
        :type coordinates: array
        :param bias: 2D array with a grid of bias
        :type bias: array
        """
        import itertools

        pnn = self._vector_eval_nonlinear_pnn(coordinates)
        bias = np.hstack([np.full((len(coordinates), 1), 1), bias])
        prod = np.array(list(itertools.combinations_with_replacement(np.arange(bias.shape[1]), r=2)))

        pgal_auto = np.zeros((pnn.shape[1], pnn.shape[2]))
        pgal_cross = np.zeros((pnn.shape[1], pnn.shape[2]))
        for i in range(len(pnn)):
            fac = 2 if prod[i, 0] != prod[i, 1] else 1
            pgal_auto += (bias[:, prod[i, 0]] * bias[:, prod[i, 1]] * fac * pnn[i, :, :].T).T
            if i < 5:
                pgal_cross += (bias[:, i] * pnn[i, :, :].T).T

        return self.nonlinear_emulator['k'], pgal_auto, pgal_cross

    def Hubble(self, omega_cold=None, omega_matter=None, omega_baryon=None,
               sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
               w0=None, wa=None, expfactor=None, k=None, **kwargs):
        """Hubble function in km / s / Mpc

        Warning: neutrino and dynamical dark energy not yet implemented
        Warning: no radiation included
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}
        O_m = kwargs['omega_matter'] if kwargs['omega_matter'] is not None else kwargs['omega_cold']
        O_Lambda = 1 - O_m
        return kwargs['hubble'] * 100 * np.sqrt(O_m / expfactor**3 + O_Lambda)

    def comoving_radial_distance(self, omega_cold=None, omega_matter=None, omega_baryon=None,
                                 sigma8_cold=None, A_s=None, hubble=None, ns=None, neutrino_mass=None,
                                 w0=None, wa=None, expfactor=None, k=None, **kwargs):
        """Comoving radial distance in Mpc/h

        Warning: neutrino and dynamical dark energy not yet implemented
        Warning: no radiation included
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key] for key in set(list(_kwargs.keys())) - set(['self'])}
        hpars = copy.deepcopy(kwargs)
        del hpars['expfactor']

        nz = int(1e4)
        z_in = 0
        z_fin = 1 / kwargs['expfactor'] - 1
        z_vector = np.linspace(z_in, z_fin, nz)
        a_vector = 1 / (1 + z_vector)
        H_vector = self.Hubble(expfactor=a_vector, **hpars)
        return 3e3 * hpars['hubble']*100 / np.trapz(H_vector, x=z_vector)


def load_lpt_emu(verbose=True):
    """Loads in memory the lpt emulator described in Aricò et al. 2021.

    :return: a dictionary containing the emulator object
    :rtype: dict
    """

    if verbose:
        print('Loading l-bias lpt emulator...')

    basefold = os.path.dirname(os.path.abspath(__file__))

    old_names = [(basefold + '/' + "lpt_emulator")]
    for old_name in old_names:
        if os.path.exists(old_name):
            import shutil
            shutil.rmtree(old_name)

    emulator_name = (basefold + '/' +
                     "lpt_emulator_v2.0.0")

    if (not os.path.exists(emulator_name)):
        import urllib.request
        import tarfile
        import ssl
        ssl._create_default_https_context = ssl._create_unverified_context
        print('Downloading emulator data (34 Mb)...')
        urllib.request.urlretrieve(
            'https://bacco.dipc.org/lpt_emulator_v2.0.0.tar',
            emulator_name + '.tar',
            MyProgressBar())
        tf = tarfile.open(emulator_name+'.tar', 'r')
        tf.extractall(path=basefold)
        tf.close()
        os.remove(emulator_name + '.tar')

    customs={"accuracy_01": accuracy_exp_01,
             "accuracy_005": accuracy_exp_005,
             "mean_absolute_exp_percentage_error":mean_absolute_exp_percentage_error}

    emulator = {}
    emulator['emu_type'] = 'nn'

    emulator['model'] = []
    emulator['sub'] = []
    emulator['scaler'] = []
    for n in range(15):
        i_emulator_name = f'{emulator_name}/lpt_emu_field{n}'

        file_to_read = open(f"{i_emulator_name}/details.pickle", "rb")
        nn_details = pickle.load(file_to_read)

        emulator['model'].append(load_model(i_emulator_name, custom_objects=customs))
        emulator['scaler'].append(nn_details['scaler'])
        emulator['sub'].append(nn_details['subtract'])

    emulator['k'] = nn_details['kk']
    emulator['keys'] = ['omega_cold', 'sigma8_cold', 'omega_baryon', 'ns',
                        'hubble', 'neutrino_mass', 'w0', 'wa', 'expfactor']
    emulator['bounds'] = nn_details['bounds']

    if verbose:
        print('L-bias lpt emulator loaded in memory.')

    return emulator

def load_nonlinear_lbias_emu(emu_type='nn', verbose=True):
    """Loads in memory the nonlinear emulator described in Zennaro et al. 2021.

    :param emu_type: type of emulator, can be 'gp' for the gaussian process, ot
                 'nn' for the neural network
    :type emu_type: str

    :return: a dictionary containing the emulator object
    :rtype: dict
    """
    if verbose:
        print('Loading non-linear l-bias emulator...')

    basefold = os.path.dirname(os.path.abspath(__file__))

    emulator_name = (basefold + '/' +
                         "lbias_emulator")

    if (not os.path.exists(emulator_name)):
        import urllib.request
        import tarfile
        import ssl
        ssl._create_default_https_context = ssl._create_unverified_context
        print('Downloading emulator data (34Mb)...')
        urllib.request.urlretrieve(
            'https://bacco.dipc.org/lbias_emulator.tar',
            emulator_name + '.tar',
            MyProgressBar())
        tf = tarfile.open(emulator_name+'.tar', 'r')
        tf.extractall(path=basefold)
        tf.close()
        os.remove(emulator_name + '.tar')
    emulator = {}
    emulator['emu_type'] = 'nn'
    emulator['model'] = []
    for n in range(15):
        i_emulator_name = f'{emulator_name}/lbias_emu_field{n}'
        emulator['model'].append(load_model(i_emulator_name))

    with open(emulator_name + '/lbias_emu.pickle', 'rb') as f:
        emulator['scaler'] = pickle.load(f)
        emulator['npca'] = pickle.load(f)
        emulator['k'] = pickle.load(f)
        _ = pickle.load(f) # components
        emulator['rotation'] = pickle.load(f)
        emulator['bounds'] = pickle.load(f)
    emulator['keys'] = ['omega_cold', 'sigma8_cold', 'omega_baryon', 'ns',
                        'hubble', 'neutrino_mass', 'w0', 'wa', 'expfactor']

    if verbose:
        print('Nonlinear l-bias emulator loaded in memory.')
    return emulator
