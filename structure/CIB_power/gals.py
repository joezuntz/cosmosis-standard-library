'''

Galaxies
Part of a series of external utils that creates kernels for Limber integrals. This is for CMB lensing.

You want to return a spline which is what is needed for limber.


'''
import numpy as np
import scipy

# ============================================================

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

# ============================================================
# ============================================================

class kern_nvss():

    def dndz_nvss_szd(z):
        '''   # Eq (21) from Smith, Zahn, Dore (2007) 0705.3980
        # Normalization to make scipy.integrate.quad( dndz_nvss_szd, 0, 1000 ) = 1 '''

    if z < 1.1:
        ret = np.exp(-(z - 1.1) ** 2 / (2 * 0.8 ** 2))
    else:
        ret = np.exp(-(z - 1.1) ** 2 / (2 * 0.3 ** 2))
    return ret * 0.82708492

    def __init__(self, omm, h0, zdndz, b=1.7):
        """ Redshift kernel for NVSS.
               cosmo - quickspec.cosmo.lcdm object describing the cosmology.
               dndz  - function dndz(z) which returns # galaxies per steradian per redshift.
               b     - linear bias parameter.
        """

        # TODO we can add shot noise if needed

        # nlgg    = 1./159000. * np.ones(lmax_dg+1) # NVSS shot noise in gal/sr, from Fig. 1 of 0705.3980

        self.b = b
    # def w_lxz(self, l, x, z):
    #     return self.b * (self.cosmo.H_z(z) * 1.e3 / units.c) * self.dndz(z)

        wa = np.zeros(np.size(zdndz))
        dndz= np.zeros(np.size(zdndz))
        zmin = zdndz[0]
        zmax = zdndz[np.size(zdndz) - 1]
        self.zmin = zmin
        self.zmax = zmax
        # TODO move this interpolation to the settings
        for i, z in enumerate(zdist):
            dndz[i] = dndz_nvss_szd(z)

        dndzfun = interp1d(zdndz, dndz[0, :])
        norm = scipy.integrate.quad(dndzfun, self.zmin, self.zmax)[0]
        wa = dndz[0, :] * b

        # Think about the H(z)/c factor which is a dchi/dz
        was = interp1d(zdndz, wa)

# ============================================================
# ============================================================


