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

        # for i, z in enumerate(zdist):
        #     wa[i] = 3. * omm / 3000. ** 2 * util.tcmb * 1.e6 * \
        #         (dzispline.derivative(1)(z) * (1. + z) / dzispline(z) + 1.)

        # # Think about the H(z)/c factor which is a dchi/dz
        # self.w_interp = scipy.interpolate.interp1d(zdist, wa)

    def w_lxz(self, l, x, z):
        return 3. * (1. + z) * self.omm / (3000. ** 2) * (x / (l+0.5)) ** 2 * \
            (self.dzispline.derivative(1)(z)) / self.dzispline(z)

# ============================================================
# ============================================================
