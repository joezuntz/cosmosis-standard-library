# Based on https://arxiv.org/abs/2011.11320 and https://arxiv.org/abs/1710.11129

from cosmosis.datablock import names
from cosmosis.gaussian_likelihood import SingleValueGaussianLikelihood

luna_mn = 0.02195
luna_sd = 0.00022

cooke18_theory_mn = 0.02166
cooke18_theory_sd = 0.00019

invcov_weighted_mn = (luna_mn/luna_sd**2 + cooke18_theory_mn/cooke18_theory_sd**2) / (1/luna_sd**2 + 1/cooke18_theory_sd**2)

sd_stat = luna_sd
sd_sys = luna_mn - invcov_weighted_mn

class BBNLikelihood(SingleValueGaussianLikelihood):
    # The mean and standard deviation of the BBN measurements.
    # The user can over-ride these in the ini file if desired
    
    mean = luna_mn
    sigma = np.sqrt(sd_stat**2 + sd_sys**2)

    # The value (either chosen by the sampler or computed
    # somewhere in the pipeline) that this Likelihood uses.
    # If we needed to operate on a derived quantity we could
    # have redefined the method block_theory_points instead
    section = names.cosmological_parameters

    # Where we should save the likelihood
    like_name = "bbn"

    def extract_theory_points(self, block):
        omega_b = block[self.section, 'omega_b']
        h = block[self.section, 'h0']
        ombh2 = omega_b * h**2
        return ombh2


setup, execute, cleanup = BBNLikelihood.build_module()
