from __future__ import print_function
from cosmosis.datablock import names
from cosmosis.gaussian_likelihood import SingleValueGaussianLikelihood

cosmo = names.cosmological_parameters
sz = "sz"


class SZXLikelihood(SingleValueGaussianLikelihood):
    # The mean and standard deviation of the Planck SZ measurement.
    # The user can over-ride these in the ini file if desired
    mean = 0.764
    sigma = 0.025
    fid_omega = 0.27

    # Where we should save the likelihood
    like_name = "sz"

    def __init__(self, options):
        super(SZXLikelihood, self).__init__(options)
        # Allow the user to override the fiducial omega
        # as well as the mean and sigma
        if options.has_value("fid_omega"):
            self.fid_omega = options["fid_omega"]
        print("Fiducial omega_m for SZ = ", self.fid_omega)

    def extract_theory_points(self, block):
        # compute the derived quantity we need for this
        # likelihood.  It is a combination of omega_m and
        # sigma_8
        sigma8 = block[cosmo, "sigma_8"]
        omega_m = block[cosmo, "omega_m"]

        # The original measurmen
        x = sigma8 * (omega_m / self.fid_omega)**0.3
        block[sz, "x"] = x
        return x


setup, execute, cleanup = SZXLikelihood.build_module()
