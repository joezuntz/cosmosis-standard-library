import os
import pickle
from cosmosis.datablock import names, SectionOptions
from hierarc.Likelihood.lens_sample_likelihood import LensSampleLikelihood
from astropy.cosmology import w0waCDM
from lenstronomy.Cosmo.cosmo_interp import CosmoInterp


class TDCOSMOlenses:
    def __init__(self, options):
        dir_path = os.path.dirname(__file__)

        self.data_sets = options.get_string("data_sets", default='tdcosmo7')
        self._num_distribution_draws = options.get_int("num_distribution_draws", default=200)
        self._distances_computation_module = options.get_string("distances_computation_module", default='astropy')

        # 7 TDCOSMO lenses
        file = open(os.path.join(dir_path, 'tdcosmo7_likelihood_processed.pkl'), 'rb')
        tdcosmo7_likelihood_processed = pickle.load(file)
        file.close()

        # 33 SLACS lenses with SDSS spectroscopy
        file = open(os.path.join(dir_path, 'slacs_sdss_likelihood_processed.pkl'), 'rb')
        slacs_sdss_likelihood_processed = pickle.load(file)
        file.close()

        # 5 SLACS with IFU
        file = open(os.path.join(dir_path, 'slacs_ifu_likelihood_processed.pkl'), 'rb')
        slacs_ifu_likelihood_processed = pickle.load(file)
        file.close()

        # here we update each individual lens likelihood configuration with the setting of the Monte-Carlo marginalization over hyper-parameter distributions
        for lens in tdcosmo7_likelihood_processed:
            lens['num_distribution_draws'] = self._num_distribution_draws
        for lens in slacs_sdss_likelihood_processed:
            lens['num_distribution_draws'] = self._num_distribution_draws
        for lens in slacs_ifu_likelihood_processed:
            lens['num_distribution_draws'] = self._num_distribution_draws

        # ====================
        # TDCOSMO 7 likelihood
        # ====================

        # hear we build a likelihood instance for the sample of 7 TDCOSMO lenses,
        lens_list = []
        if 'tdcosmo7' in self.data_sets:
            lens_list += tdcosmo7_likelihood_processed
        if 'SLACS_SDSS' in self.data_sets:
            lens_list += slacs_sdss_likelihood_processed
        if 'SLACS_IFU' in self.data_sets:
            lens_list += slacs_ifu_likelihood_processed

        assert len(
            lens_list) > 0, "Data not found ! Add at least one of those 3 data sets 'tdcosmo7', 'SLACS_SDSS' or 'SLACS_IFU'"
        # choose which likelihood you want here:
        self._likelihood = LensSampleLikelihood(lens_list)

        # choose if you want the full astropy distance calculation or a interpolated version of it (for speed-up)
        self._interpolate_distances_type = 'None'

    def cosmosis_cosmo_2_astropy_cosmo(self, block):
        """

        :param cosmosis_cosmo: cosmosis cosmology object
        :return ~astropy.cosmology equivalent cosmology object
        """
        H0 = block['cosmological_parameters', 'h0'] * 100
        om = block['cosmological_parameters', 'omega_m']
        ok = block['cosmological_parameters', 'omega_k']
        ob = block['cosmological_parameters', 'omega_b']
        w0 = block['cosmological_parameters', 'w']
        wa = block['cosmological_parameters', 'wa']
        mnu = block['cosmological_parameters', 'mnu']
        nnu = block['cosmological_parameters', 'nnu']

        ol = 1 - om - ok

        if self._distances_computation_module == 'astropy':
            # we are using standard astropy cosmology for distance compuation
            cosmo = w0waCDM(H0=H0, Om0=om, Ode0=ol, Ob0=ob, w0=w0, wa=wa, m_nu=mnu, Neff=nnu)
        elif self._distances_computation_module == 'CosmoInterp':
            # we are using an interpolated version of the standard astropy cosmology (for speed-up)
            cosmo = w0waCDM(H0=H0, Om0=om, Ode0=ol, Ob0=ob, w0=w0, wa=wa, m_nu=mnu, Neff=nnu)
            cosmo = CosmoInterp(cosmo=cosmo, z_stop=5, num_interp=100)
        elif self._distances_computation_module == 'camb':

            # we use the camb distances
            z_bg = block['distances', 'z']
            D_A = block['distances', 'D_A']
            # todo: compute K = Ok * c^2/H0^2 as dimensionless units
            cosmo = CosmoInterp(ang_dist_list=D_A, z_list=z_bg, Ok0=ok, K=None,)
        else:
            raise ValueError()

        return cosmo

    def likelihood(self, block):
        cosmo = self.cosmosis_cosmo_2_astropy_cosmo(block)

        # here the additional parameters required to evaluate the likelihood in accordance with TDCOSMO IV Table 3
        lambda_mst = block['nuisance_strong_lensing', 'lambda_mst']
        log_lambda_mst_sigma = block['nuisance_strong_lensing', 'log_lambda_mst_sigma']
        lambda_mst_sigma = 10**log_lambda_mst_sigma
        
        alpha_lambda = block['nuisance_strong_lensing', 'alpha_lambda']

        # We will define this parameter in the block in log space because the prior is uniform in log_ space. 
        # a_ani = block['nuisance_strong_lensing', 'a_ani']
        # a_ani_sigma = block['nuisance_strong_lensing', 'a_ani_sigma']
        
        log_a_ani = block['nuisance_strong_lensing', 'log_a_ani']
        a_ani = 10**log_a_ani
        
        log_a_ani_sigma = block['nuisance_strong_lensing', 'log_a_ani_sigma']
        a_ani_sigma = 10**log_a_ani_sigma


        kwargs_lens_test = {'lambda_mst': lambda_mst,  # mean in the internal MST distribution
                            'lambda_mst_sigma': lambda_mst_sigma,  # Gaussian sigma of the distribution of lambda_mst
                            'alpha_lambda': alpha_lambda,  # slope of lambda_mst with r_eff/theta_E
                            }
        kwargs_kin_test = {'a_ani': a_ani,  # mean a_ani anisotropy parameter in the OM model
                           'a_ani_sigma': a_ani_sigma,  # sigma(a_ani)⟨a_ani⟩ is the 1-sigma Gaussian scatter in a_ani
                           }

        logl = self._likelihood.log_likelihood(cosmo=cosmo, kwargs_lens=kwargs_lens_test, kwargs_kin=kwargs_kin_test)

        return float(logl)


def setup(options):
    options = SectionOptions(options)
    return TDCOSMOlenses(options)


def execute(block, config):
    like = config.likelihood(block)
    block[names.likelihoods, "TDCOSMO_like"] = like
    return 0

# todo : Here implement a custom class that have the same method as astropy cosmology, but overwrite the angular diameter computation
# class CustomCosmo(w0waCDM):
#     def __init__(z_bkg, D_A):
#         self.z_bkg = z_bkg
#         self.D_A = D_A
#         self.angular_distance = interp1d(self.z_bkg, self.D_A, kind='linear')
#     def Da(z):
#           return self.angular_distance(z)
#     def calculate_Dds(self, ok, K, zd, zs, Dd, Ds):
#
#         # to compute Dds from Dd and Ds, first need to figure out whether the universe is flat
#
#         if (np.fabs(ok) < 1.e-6):  # flat Universe, so can use simple formula
#             Dds = ((1. + zs) * Ds - (1 + zd) * Dd) / (1. + zs)
#
#         elif (K > 0):
#             chis = np.arcsin(Ds * (1. + zs) * np.sqrt(K)) / np.sqrt(K)
#             chid = np.arcsin(Dd * (1. + zd) * np.sqrt(K)) / np.sqrt(K)
#             chids = chis - chid
#             Dds = (1. / (1 + zs)) * (1. / np.sqrt(K)) * np.sin(np.sqrt(K) * chids)
#
#         else:  # K<0
#             chis = np.arcsinh(Ds * (1. + zs) * np.sqrt(-K)) / np.sqrt(-K)
#             chid = np.arcsinh(Dd * (1. + zd) * np.sqrt(-K)) / np.sqrt(-K)
#             chids = chis - chid
#             Dds = (1. / (1 + zs)) * (1. / np.sqrt(-K)) * np.sinh(np.sqrt(-K) * chids)
#
#         return Dds
