import os
import pickle
from cosmosis.datablock import names, SectionOptions
from hierarc.Likelihood.lens_sample_likelihood import LensSampleLikelihood
from astropy.cosmology import LambdaCDM
from lenstronomy.Cosmo.cosmo_interp import CosmoInterp
from scipy.interpolate import interp1d

class TDCOSMOlenses:
    def __init__(self, options):
        dir_path = os.path.dirname(__file__)

        #todo: unpack options to select the cosmology. Default is LCDM
        self.cosmo_model = options.get_string("cosmo_model", default='LambdaCDM')
        if self.cosmo_model == 'LambdaCDM':
            self.cosmology = LambdaCDM
        else:
            raise NotImplementedError()

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


        num_distribution_draws = 200  # number of draws from the hyper-parameter distribution in computing the Monte Carlo integral marginalization

        # here we update each individual lens likelihood configuration with the setting of the Monte-Carlo marginalization over hyper-parameter distributions
        for lens in tdcosmo7_likelihood_processed:
            lens['num_distribution_draws'] = num_distribution_draws
        for lens in slacs_sdss_likelihood_processed:
            lens['num_distribution_draws'] = num_distribution_draws
        for lens in slacs_ifu_likelihood_processed:
            lens['num_distribution_draws'] = num_distribution_draws

        # ====================
        # TDCOSMO 7 likelihood
        # ====================

        # hear we build a likelihood instance for the sample of 7 TDCOSMO lenses,
        self.data_sets = options.get_string("data_sets", default='tdcosmo7')
        lens_list = []
        if 'tdcosmo7' in self.data_sets:
            lens_list += tdcosmo7_likelihood_processed
        if 'SLACS_SDSS' in self.data_sets:
            lens_list += slacs_sdss_likelihood_processed
        if 'SLACS_IFU' in self.data_sets:
            lens_list += slacs_ifu_likelihood_processed

        assert len(lens_list)>0, "Data not found ! Add at least one of those 3 data sets 'tdcosmo7', 'SLACS_SDSS' or 'SLACS_IFU'"
        #choose which likelihood you want here:
        self.likelihood = LensSampleLikelihood(lens_list)

        #choose if you want the full astropy distance calculation or a interpolated version of it (for speed-up)
        self._interpolate_distances_type = 'None'

    def cosmosis_cosmo_2_astropy_cosmo(self, block):
        """

        :param cosmosis_cosmo: cosmosis cosmology object
        :return ~astropy.cosmology equivalent cosmology object
        """
        h = block['cosmological_parameters', 'h0']
        om = block['cosmological_parameters', 'omega_m']
        ol = block['cosmological_parameters', 'omega_l'] #todo: check this, not sure how parameters are named

        if self._interpolate_distances_type == 'CosmoInterp':
            cosmo = self.cosmology(H0=h, Om0=om, Ode0=ol)
            cosmo = CosmoInterp(cosmo=cosmo, z_stop=5, num_interp=100)
        elif self._interpolate_distances_type == 'None':
            #we are falling back to astropy cosmology (might be slow)
            cosmo = self.cosmology(H0=h, Om0=om, Ode0=ol)
        elif self._interpolate_distances_type == 'cosmosis':
            #use the distance redshift relation of cosmosis
            z_bg = block['distances', 'z']
            D_A = block['distances', 'D_A']
            angular_distance = interp1d(z_bg, D_A, kind='linear')
            # todo to be implemented as in the holicow module
            raise NotImplementedError()
        else:
            raise ValueError()



        return cosmo


    def likelihood(self, block):
        cosmo = self.cosmosis_cosmo_2_astropy_cosmo(block)

        # here the additional parameters required to evaluate the likelihood in accordance with TDCOSMO IV Table 3
        # todo: add those as nuisance parameters, they are fixed for now
        kwargs_lens_test = {'lambda_mst': 1.,  # mean in the internal MST distribution
                        'lambda_mst_sigma': 0.04,  # Gaussian sigma of the distribution of lambda_mst
                        'alpha_lambda': 0,  # slope of lambda_mst with r_eff/theta_E
                       }
        kwargs_kin_test = {'a_ani': 1.5,  # mean a_ani anisotropy parameter in the OM model
                       'a_ani_sigma': 0.3,  # sigma(a_ani)⟨a_ani⟩ is the 1-sigma Gaussian scatter in a_ani
                      }

        logl = self.likelihood.log_likelihood(cosmo=cosmo, kwargs_lens=kwargs_lens_test, kwargs_kin=kwargs_kin_test)

        return float(logl)

def setup(options):
    options = SectionOptions(options)
    return TDCOSMOlenses(options)

def execute(block, config):
    like = config.likelihood(block)
    block[names.likelihoods, "TDCOSMO_like"] = like
    return 0

