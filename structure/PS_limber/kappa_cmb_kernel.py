'''

KAPPA CMB

Part of a series of external utils that creates kernels for Limber integrals. This one is for CMB lensing.

You want to return a spline function W(l,chi,z) with l multipole chi comiving distance z redsfhit which is what is needed for limber.

EVERYTHING IS IN h UNITS

'''

import numpy as np
import scipy


def chiint(z, omegam, h0):
    """
    Comoving distance integral in case CAMB did not give you the evolution till high redshifts
    """
    # 3000 is the c in H/c in h units
    chiint = 3000. / np.sqrt(omegam * (1. + z) ** 3 + (1. - omegam))
    return chiint


class kern():

            '''
        KAPPA CMB KERNEL (h units):

        Args:

            zdist: redshift distribution of the spline
            omm: Omega matter
            h0: Hubble constant
            xlss: Last scattering surface comoving distance

        Return:

            kern().w_lxz: kernel for limber integral


       '''

    def __init__(self, zdist, omm, h0, xlss):
        wb = np.zeros(np.size(zdist))
        # use the z's from the P(k,z) array
        zmax = zdist[np.size(zdist) - 1]
        zmin = zdist[0]
        zmax = zdist[np.size(zdist) - 1]
        self.h0 = h0
        self.omm = omm
        self.zmin = zmin
        self.zmax = zmax
        self.xlss = xlss
        chicmb = xlss

    def w_lxz(self, l, x, z):
        '''
        KAPPA CMB KERNEL (h units):

        wcmb = 1.5*omegam*h0**2*(1.+z)*chi*(1-chi(z)/chi(zcmb))
        '''
        chiz = scipy.integrate.quad(chiint, 0., z, args=(self.omm, self.h0))[0]
        if (z < self.zmax / 1.0001):
            tmp = (1. - chiz / self.xlss)

        return 1.5 * self.omm * (1. + z) * chiz * tmp / 3000. ** 2
