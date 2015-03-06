import numpy as np

import scipy.integrate
import scipy.interpolate

import quickspec as qs

#--

class halo_model(qs.mps.mps):
    """ base class for encapsulating a halo model for large scale structure (see Cooray and Sheth 2002 http://arxiv.org/abs/astro-ph/0206508) """
    def __init__(self, mass_function, halo_profile, hod, p_lin):
        self.mass_function = mass_function
        self.halo_profile  = halo_profile
        self.hod           = hod
        self.p_lin         = p_lin

        self.cosmo         = p_lin.cosmo

    def p_kz(self, k, z):
        return self.p_kz_1h(k, z) + self.p_kz_2h(k, z)

    def p_kz_1h(self, k, z):
        def integrand(logm):
            m = np.exp(logm)
            return m * self.mass_function.dndM_mz(m,z) * self.hod.hod_1h(m, z) * self.halo_profile.u_km(k,m,z)**2
        return scipy.integrate.quad( integrand, np.log(self.mass_function.Mmin), np.log(self.mass_function.Mmax) )[0]

    def p_kz_2h(self, k, z):
        def integrand(logm):
            m = np.exp(logm)
            return m * self.mass_function.dndM_mz(m,z) * self.hod.hod_2h(m, z) * self.mass_function.b_mz(m, z) * self.halo_profile.u_km(k,m,z)
        return scipy.integrate.quad( integrand, np.log(self.mass_function.Mmin), np.log(self.mass_function.Mmax) )[0]**2 * self.p_lin.p_kz(k, z)

class halo_model_cache(halo_model):
    """ subclass of halo_model which precomputes the 1h and 2h power spectrum terms over a grid of wavenumbers k and redshifts z, for interpolation later. """
    def __init__(self, mod,
                 kmin=1.e-3, kmax=1., npts=10,
                 zs = np.array([0.0, 0.5, 1., 1.5, 2., 3., 4., 6., 8., 12., 16., 32., 64., 128., 256., 512., 1300.])):

        self.cosmo = mod.p_lin.cosmo

        self.model = mod
        self.kmin  = kmin
        self.kmax  = kmax

        self.vec_z   = zs
        self.vec_k   = np.logspace( np.log10(kmin), np.log10(kmax), (np.log10(kmax) - np.log10(kmin))*npts )
        self.vec_lnk = np.log(self.vec_k)

        self.mat_lnp_kz_1h = np.zeros( (len(self.vec_k), len(self.vec_z)) )
        self.mat_lnp_kz_2h = np.zeros( (len(self.vec_k), len(self.vec_z)) )

        for ik, k in qs.util.enumerate_progress(self.vec_k, label="halo::model_cache::init"):
            for iz, z in enumerate(self.vec_z):
                self.mat_lnp_kz_1h[ik, iz] = np.log(mod.p_kz_1h(k, z))
                self.mat_lnp_kz_2h[ik, iz] = np.log(mod.p_kz_2h(k, z) / mod.p_lin.p_kz(k, z))
        self.spl_lnp_kz_1h = scipy.interpolate.RectBivariateSpline( self.vec_lnk, self.vec_z, self.mat_lnp_kz_1h, kx=3, ky=3, s=0)
        self.spl_lnp_kz_2h = scipy.interpolate.RectBivariateSpline( self.vec_lnk, self.vec_z, self.mat_lnp_kz_2h, kx=3, ky=3, s=0)

        self.kmin = np.min(self.vec_k)
        self.kmax = np.max(self.vec_k)
        self.zmin = np.min(self.vec_z)
        self.zmax = np.max(self.vec_z)

    def p_kz_1h(self, k, z):
        assert( np.all(k >= self.kmin) )
        assert( np.all(k <= self.kmax) )
        assert( np.all(z >= self.zmin) )
        assert( np.all(z <= self.zmax) )

        return np.exp( self.spl_lnp_kz_1h.ev(np.log(k), z) )

    def p_kz_2h(self, k, z):
        assert( np.all(k >= self.kmin) )
        assert( np.all(k <= self.kmax) )
        assert( np.all(z >= self.zmin) )
        assert( np.all(z <= self.zmax) )

        return np.exp( self.spl_lnp_kz_2h.ev(np.log(k), z) ) * self.model.p_lin.p_kz(k, z)

# --

model = halo_model
model_cache = halo_model_cache
