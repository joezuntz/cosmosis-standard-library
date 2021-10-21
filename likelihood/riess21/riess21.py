from cosmosis.datablock import names
from cosmosis.gaussian_likelihood import SingleValueGaussianLikelihood


class Riess21Likelihood(SingleValueGaussianLikelihood):
    # The mean and standard deviation of the Riess11 measurements.
    # The user can over-ride these in the ini file if desired
    mean = 0.732
    sigma = 0.013

    # The value (either chosen by the sampler or computed
    # somewhere in the pipeline) that this Likelihood uses.
    # If we needed to operate on a derived quantity we could
    # have redefined the method block_theory_points instead
    section = names.cosmological_parameters
    name = "h0"

    # Where we should save the likelihood
    like_name = "riess21"


setup, execute, cleanup = Riess21Likelihood.build_module()
