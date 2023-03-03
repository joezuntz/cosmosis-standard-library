"""
The likelihood module for Pantheon+ SNe combined with SH0ES Cepheids

"""

from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
import os
import numpy as np
import pandas as pd
import gzip
import pantheon_covariance_io as pan_io

# Default is to use SN data from https://arxiv.org/abs/2202.04077
# and Cepheids data from https://arxiv.org/abs/2112.04510
default_data_file = os.path.join(os.path.split(__file__)[0], "Pantheon+SH0ES.dat")
default_covmat_file = os.path.join(os.path.split(__file__)[0], "Pantheon+SH0ES_STAT+SYS.cov_compressed.gz")


class PantheonLikelihood(GaussianLikelihood):
    x_section = names.distances
    x_name = "z"
    y_section = names.distances
    y_name = "D_A"
    like_name = "pantheon"

    def build_data(self):
        """
        Run once at the start to load in the data vectors.

        Returns x, y where x is the independent variable (redshift in this case)
        and y is the Gaussian-distribured measured variable (magnitude in this case).

        """
        filename = self.options.get_string("data_file", default=default_data_file)

        # Decide whether or not to include the extra SHOES data
        self.include_shoes = self.options.get_bool("include_shoes")
        if self.include_shoes:
            print("Using PantheonPlus and SHOES")
        else:
            print("Using PantheonPlus only (not include SHOES)")

        print("Loading data from {}".format(filename))
        data = pd.read_csv(filename,delim_whitespace=True)
        self.origlen = len(data)
        
        if self.include_shoes:
            self.ww = (data['zHD']>0.01) | (np.array(data['IS_CALIBRATOR'],dtype=bool))
            self.is_calibrator = data['IS_CALIBRATOR'][self.ww].astype(bool)
            self.cepheid_distance = data['CEPH_DIST'][self.ww]
        else:
            self.ww = (data['zHD']>0.01)


        #use the vpec corrected redshift for zCMB 
        self.zCMB = data['zHD'][self.ww] 
        self.zHEL = data['zHEL'][self.ww]
        self.m_obs = data['m_b_corr'][self.ww]

        return self.zCMB, self.m_obs


    def build_covariance(self):
        """Run once at the start to build the covariance matrix for the data"""
        filename = self.options.get_string("covmat_file", default=default_covmat_file)
        print("Loading covariance from {}".format(filename))

        # Read the covariance
        if filename.endswith('.gz'):
            C = pan_io.read_compressed(filename)
        else:
            C = pan_io.read_uncompressed(filename)

        C = C[self.ww][:, self.ww]

        # Return the covariance; the parent class knows to invert this
        # later to get the precision matrix that we need for the likelihood.
        return C

    def extract_theory_points(self, block):
        """
        Run once per parameter set to extract the mean vector that our
        data points are compared to.  For the Hubble flow set, we compare to the 
        cosmological model. For the calibrators we compare to the Cepheid distances.
        """
        import scipy.interpolate

        # Pull out theory mu and z from the block.
        theory_x = block[self.x_section, self.x_name]
        theory_y = block[self.y_section, self.y_name]
        theory_ynew = self.zCMB * np.nan

        # Interpolation function of theory so we can evaluate at redshifts of the data
        f = scipy.interpolate.interp1d(theory_x, theory_y, kind=self.kind)
        
        if self.include_shoes:
            #Here we use the Cepheid host distances as the "theory". 
            #This is the 2nd rung of the distance ladder and is what calibrates M
            theory_ynew[self.is_calibrator] = self.cepheid_distance[self.is_calibrator]

            # Actually do the interpolation at the data redshifts
            zcmb = self.zCMB[~self.is_calibrator]
            zhel = self.zHEL[~self.is_calibrator]
            theory_ynew[~self.is_calibrator] = 5.0*np.log10((1.0+zcmb)*(1.0+zhel)*np.atleast_1d(f(zcmb)))+25.
        else:
            zcmb = self.zCMB
            zhel = self.zHEL
            theory_ynew = 5.0*np.log10((1.0+zcmb)*(1.0+zhel)*np.atleast_1d(f(zcmb)))+25.

        # Add the absolute supernova magnitude and return
        M = block[names.supernova_params, "M"]
        return theory_ynew + M


# This takes our class and turns it into 
setup, execute, cleanup = PantheonLikelihood.build_module()
