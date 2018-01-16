from cosmosis.datablock import names
from cosmosis.gaussian_likelihood import SingleValueGaussianLikelihood


class Riess11Likelihood(SingleValueGaussianLikelihood):
    # The mean and standard deviation of the Riess11 measurements.
    # The user can over-ride these in the ini file if desired
    mean = 0.7324
    sigma = 0.0174

    # The value (either chosen by the sampler or computed
    # somewhere in the pipeline) that this Likelihood uses.
    # If we needed to operate on a derived quantity we could
    # have redefined the method block_theory_points instead
    section = names.cosmological_parameters
    name = "h0"

    # Where we should save the likelihood
    like_name = "riess16"


setup, execute, cleanup = Riess11Likelihood.build_module()
