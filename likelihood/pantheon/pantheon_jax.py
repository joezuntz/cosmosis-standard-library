"""
The likelihood module for the Pantheon dataset discards
much of the JLA/CosmoMC SN machinery, which is not needed
for Pantheon because the systematic errors in that dataset
have been subsumed into a single systematic covariance matrix.

In consequence almost all the terms in those codes are zero.

"""

# We don't need cosmosis to jaxify our module as it is already a jax module
cosmosis_jax = False

from cosmosis.gaussian_likelihood_jax import GaussianLikelihood
from cosmosis.datablock import names
import os
import numpy as np
import sys
import jax

# Default is to use the binned version of the data since it's much faster
# You can also downloaded and run the full data if you like, and set the data_file
# and covmat_file parameters in the ini file.
default_data_file = os.path.join(os.path.split(__file__)[0], "lcparam_DS17f.txt")
default_covmat_file = os.path.join(os.path.split(__file__)[0], "sys_DS17f.txt")


class PantheonLikelihood(GaussianLikelihood):
    x_section = names.distances
    x_name = "z"
    y_section = names.distances
    y_name = "mu"
    like_name = "pantheon"


    def build_data(self):
        """
        Run once at the start to load in the data vectors.

        Returns x, y where x is the independent variable (redshift in this case)
        and y is the Gaussian-distribured measured variable (magnitude in this case).

        """
        from mpi4py.MPI import COMM_WORLD as comm

        if comm.rank == 0:
            filename = self.options.get_string("data_file", default=default_data_file)
            print("Loading Pantheon data from {}".format(filename))

            # The Pantheon data format is mostly zeros.
            # The only columns that we actually need here are the redshift,
            # magnitude, and magnitude error.
            data = np.genfromtxt(filename).T
            z = data[1]
            m_obs = data[4]

            # We will need mag_obs_err later, when building the covariance,
            # so save it for now.
            self.mag_obs_err = data[5]

            # Return this to the parent class, which will use it
            # when working out the likelihood
            print("Found {} Pantheon supernovae (or bins if you used the binned data file)".format(len(z)))

            if comm.size > 1:
                comm.bcast(z, root=0)
                comm.bcast(m_obs, root=0)
                comm.bcast(self.mag_obs_err, root=0)
        
        else:
            z = comm.bcast(None, root=0)
            m_obs = comm.bcast(None, root=0)
            self.mag_obs_err = comm.bcast(None, root=0)
            sys.stderr.write(f"Rank {comm.rank} got {len(z)} supernovae bcast\n")
        print("Pantheon redshift range: {} to {}".format(z.min(), z.max()))
        return z, m_obs

    def build_covariance(self):
        """Run once at the start to build the covariance matrix for the data"""
        filename = self.options.get_string("covmat_file", default=default_covmat_file)
        print("Loading Pantheon covariance from {}".format(filename))
        # The file format for the covariance has the first line as an integer
        # indicating the number of covariance elements, and the the subsequent
        # lines being the elements.
        # This data file is just the systematic component of the covariance - 
        # we also need to add in the statistical error on the magnitudes
        # that we loaded earlier
        f = open(filename)
        line = f.readline()
        n = int(line)
        C = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                C[i,j] = float(f.readline())

        # Now add in the statistical error to the diagonal
        for i in range(n):
            C[i,i] += self.mag_obs_err[i]**2
        f.close()

        # Return the covariance; the parent class knows to invert this
        # later to get the precision matrix that we need for the likelihood.
        return C

    def extract_theory_samples(self, block):
        "Extract relevant theory from block at sample points"
        z = block[self.x_section, self.x_name]
        mu = block[self.y_section, self.y_name]

        if z[0] == 0:
            dz = z[1] - z[0]
            dmu = mu[2] - mu[1]
            mu[0] = mu[1] + dmu/dz * (-dz)
        
        M = block[names.supernova_params, "M"]

        # This is kind of a hack.  Actually we would like to do the likelihood
        # with respect to the M parameter, but this has the same effect.
        block.put_derivative("distances", "mu", names.supernova_params, "M", np.ones_like(mu))
        return z, mu + M




# This takes our class and turns it into 
setup, execute, cleanup = PantheonLikelihood.build_module()
