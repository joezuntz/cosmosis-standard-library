'''

Galaxies
Part of a series of external utils that creates kernels for Limber integrals. This is for CMB lensing.

You want to return a spline which is what is needed for limber.


'''


class kern():

    def __init__(self, omm, h0, zdndz, dndz, b=1.):
        """ Redshift kernel for a galaxy survey.
               cosmo - quickspec.cosmo.lcdm object describing the cosmology.
               dndz  - function dndz(z) which returns # galaxies per steradian per redshift.
               b     - linear bias parameter.
        """

        self.b = b
    # def w_lxz(self, l, x, z):
    #     return self.b * (self.cosmo.H_z(z) * 1.e3 / units.c) * self.dndz(z)

        wa = np.zeros(np.size(zdndz))
        zmin = zdndz[0]
        zmax = zdndz[np.size(zdndz) - 1]
        self.zmin = zmin
        self.zmax = zmax
        # TODO move this interpolation to the settings
        dndzfun = interp1d(zdndz, dndz[0, :])
        norm = scipy.integrate.quad(dndzfun, self.zmin, self.zmax)[0]
        wa = dndz[0, :] * b
        # Think about the H(z)/c factor which is a dchi/dz
        was = interp1d(zdndz, wa)


import numpy as np
import scipy


def chiint(z, omegam, h0):
    chiint = 3000. / np.sqrt(omegam * (1. + z) ** 3 + (1. - omegam))
    return chiint


class kern():
