'''

Galaxies
Part of a series of external utils that creates kernels for Limber integrals. This is for CMB lensing.

You want to return a spline which is what is needed for limber.

A nice review of redshift models for NVSS can be find in

- 1312.0530


'''
import numpy as np
import scipy

# AUXILIARY FUNCTIONS

# NVSS DNDZ

def dndz_nvss_szd(z):
    '''   # Eq (21) from Smith, Zahn, Dore (2007) 0705.3980
    # Normalization to make scipy.integrate.quad( dndz_nvss_szd, 0, 1000 ) = 1 '''

    if z < 1.1:
        ret = np.exp(-(z - 1.1) ** 2 / (2 * 0.8 ** 2))
    else:
        ret = np.exp(-(z - 1.1) ** 2 / (2 * 0.3 ** 2))
    return ret * 0.82708492


def dndz_nvss_dezotti(z):
    '''
    de Zotti et al. (2010):
    dn/dz= 1.29 + 32.37*z - 32.89*z**2 + 11.13*z**3 - 1.25*z**4
    '''
    return 1.29 + 32.37 * z - 32.89 * z ** 2 + 11.13 * z ** 3 - 1.25 * z ** 4


def dndz_nvss_general(z, z0=0.53, alpha=0.81):
    '''
    Gamma function is a standard way to define redshift distribution of galaxies.
    1312.0530 find for NVSS

    alpha=0.81
    z0=0.53

    '''

    return (z / z0) ** alpha * np.exp(-alpha * z / z0)


# ============================================================

class kern():

    def __init__(self, omm, h0, hspline, z, dndz=dndz_nvss_szd, b=1.7):
        '''
        Galaxies KERNEL (h units):

        Args:

            zdist: redshift distribution of the spline
            dndzfun: galaxies redshift distribution
            omm: Omega matter
            h0:hubble constant
            b: Galaxies bias


        Return:

            kern().w_lxz: kernel for limber integral


       '''

        # TODO we can add shot noise if needed

        # nlgg    = 1./159000. * np.ones(lmax_dg+1) # NVSS shot noise in gal/sr, from Fig. 1 of 0705.3980

        self.b = b
    # def w_lxz(self, l, x, z):
    #     return self.b * (self.cosmo.H_z(z) * 1.e3 / units.c) * self.dndz(z)

        zmin = z[0]
        zmax = z[np.size(z) - 1]
        self.zmin = zmin
        self.zmax = zmax
        self.norm = scipy.integrate.quad(dndz, self.zmin, self.zmax)[0]
        self.dndz = dndz
        self.hspline = hspline
        # Think about the H(z)/c factor which is a dchi/dz

# ============================================================
# ============================================================

    def w_lxz(self, l, x, z):
        '''
        Galaxies KERNEL (h units):

        Notes:

        delta_g = bias_g \int dchi dn/d\chi delta_m(\chi,t)
        see eq 10 of 0705.3980 so
        w = H(z) dN/dz * b(z)

        H(z) compensate the dN/dz

        with
        '''

        return self.hspline(z) * self.dndz(z) / self.norm * self.b
