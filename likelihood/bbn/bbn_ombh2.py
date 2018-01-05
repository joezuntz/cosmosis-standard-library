from cosmosis.datablock import names
from cosmosis.gaussian_likelihood import SingleValueGaussianLikelihood


class BBNLikelihood(SingleValueGaussianLikelihood):
    # The mean and standard deviation of the BBN measurements.
    # The user can over-ride these in the ini file if desired
    mean = 0.023
    sigma = 0.002

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
