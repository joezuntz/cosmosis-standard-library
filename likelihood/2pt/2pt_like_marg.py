from __future__ import print_function
from builtins import zip
from builtins import range
from builtins import object
from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
from twopoint_cosmosis import theory_names, type_table
from astropy.io import fits
from scipy.interpolate import interp1d
import numpy as np
import twopoint
import gaussian_covariance
import os
twopt_like = __import__("2pt_like")

default_array = np.repeat(-1.0, 99)

class TwoPointPointMassLikelihood(twopt_like.TwoPointLikelihood):
        like_name='2pt'
        def __init__(self, options):
            super(TwoPointPointMassLikelihood, self).__init__(options)
        
        def build_data(self):
                #build the data as in the usual TwoPointLikelihood
                _,datavector = super(TwoPointPointMassLikelihood, self).build_data()

                #Now if requested do the point mass mode projection
                #We are essentially adding sigma_a^2 * (1/theta^2) to the covariance.
                #But this can make the inversion numerically unstable. So instead, 
                #use the beautiful Sherman-Morrison formula to amend the inverse covariance
                self.inv_cov = np.linalg.inv( self.two_point_data.covmat )
                self.sigma_a = self.options.get_double( "sigma_a", -1. )
                template_vector = np.zeros( self.inv_cov.shape[0] )
                if self.sigma_a > 0.:
                        if "gammat" in self.used_names:
                                gammat_start_ind=0
                                for s in self.two_point_data.spectra:
                                        if s.name!="gammat":
                                                gammat_start_ind+=len(s.value)
                                        else:
                                                gammat_spec=s
                                                break

                                #Now loop through lens-source pairs doing the necessary covariance doctoring....
                                for bin_pair in gammat_spec.get_bin_pairs():
                                        lens_bin,source_bin = bin_pair
                                        use = (s.bin1==lens_bin)*(s.bin2==source_bin)
                                        spec_inds, = np.where(use)
                                        dv_inds = spec_inds + gammat_start_ind
                                        theta_vals = s.angle[use]
                                        theta_m2 = theta_vals**-2
                                        template_vector[dv_inds] = theta_m2
                template_vector *= self.sigma_a
                #Now use SM to amend the inverse covariance
                uvT = np.outer(template_vector, template_vector)
                num = np.matmul( self.inv_cov, np.matmul( uvT, self.inv_cov ) )
                denom = 1 + np.dot( template_vector, np.dot( self.inv_cov, template_vector ) )
                sub_inv_cov = num/denom
                print('subtracting from inverse covariance:')
                print(sub_inv_cov)
                self.inv_cov -= sub_inv_cov
                return None, datavector

        def build_inverse_covariance(self):
                """
                Override the build_inverse_covariance method to change
                how the inverse is generated from the covariance.

                When the covariance is generated from a suite of simulations,
                for example, the simple inverse is not the best estimate.

                """
                return self.inv_cov

setup, execute, cleanup = TwoPointPointMassLikelihood.build_module()
