'''

ISW

Part of a series of external utils that creates kernels for Limber integrals. This is for CMB lensing.

You want to return a spline which is what is needed for limber.


'''
import numpy as np
import scipy
import util

# ============================================================


class kern():

    def __init__(self, h0, omm, zdist, chispline, dzispline):
        """ Redshift kernel for a galaxy survey.
               cosmo - quickspec.cosmo.lcdm object describing the cosmology.
               dndz  - function dndz(z) which returns # galaxies per steradian per redshift.
               b     - linear bias parameter.
        """

        wa = np.zeros(np.size(zdist))
        zmin = zdist[0]
        zmax = zdist[np.size(zdist) - 1]
        self.h0 = h0
        self.omm = omm
        self.zmin = zmin
        self.zmax = zmax
        self.dzispline = dzispline



    def w_lxz(self, l, x, z):

        print l,x,z,(self.dzispline.derivative(1)(z)) / self.dzispline(z)
        # yhe 1/ho**3 is coming as in CIB_kernel hall.py to compensate for the h factor.
        # TODO recheck and be sure or, formulate the code to be h independent (no division by h and check hspline the is maybe alrady divided)
        return 3. * (1. + z) * self.omm / (3000. ** 2) * (x / (l+0.5)) ** 2 * \
            (self.dzispline.derivative(1)(z)) / self.dzispline(z)/ self.h0 ** 3



# ============================================================
# ============================================================
