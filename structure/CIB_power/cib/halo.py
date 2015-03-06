import numpy as np

import scipy.integrate
import scipy.interpolate

import quickspec as qs

class hod_cib_pep():
    # HOD prescription from the Planck early CIB paper.
    # arxiv:1101.2028
    def __init__(self, mf, Mmin, asat):
        self.mf = mf

        self.Mmin  = Mmin
        self.asat  = asat
        self.Msat  = 3.3 * Mmin
        self.slogM = 0.65

        self.zvec    = np.linspace(0, 1300., 10000)
        self.lnnvec  = np.zeros( len(self.zvec) )

        for iz, z in qs.util.enumerate_progress(self.zvec, "halo::hod::sat::init"):
            tn = scipy.integrate.quad( lambda lnm : ( self.mf.dndM_mz(np.exp(lnm), z) *
                                                      self.ngal(np.exp(lnm), z) *
                                                      np.exp(lnm) ),
                                       np.log(mf.Mmin), np.log(mf.Mmax) )[0]

            self.lnnvec[iz] = np.log( max(tn, 1.3e-87) )

        self.spl_lnn = scipy.interpolate.UnivariateSpline( self.zvec, self.lnnvec, k=3, s=0 )

        self.zmin = np.min(self.zvec)
        self.zmax = np.max(self.zvec)

    def nbar(self, z):
        assert( np.all( z >= self.zmin ) )
        assert( np.all( z <= self.zmax ) )

        return np.exp( self.spl_lnn(z) )

    def ngal(self, m, z):
        return self.ncen(m, z) + self.nsat(m, z)

    def ncen(self, m, z):
        return 0.5 * ( 1.+scipy.special.erf(np.log10(m/self.Mmin)/self.slogM) )

    def nsat(self, m, z):
        return 0.5 * ( 1.+scipy.special.erf(np.log10(0.5*m/self.Mmin)/self.slogM) )*(m/self.Msat)**self.asat

    def hod_1h(self, m, z):
        nsat = self.nsat(m, z)
        ncen = self.ncen(m, z)
        nbar = self.nbar(z)
        return (2*nsat*ncen + nsat**2) / nbar**2

    def hod_2h(self, m, z):
        ngal = self.ngal(m, z)
        nbar = self.nbar(z)
        return ngal / nbar

class model_cib_x_phi(qs.mps.mps):
    def __init__(self, mass_function, halo_profile, hod, p_lin):
        self.mass_function = mass_function
        self.halo_profile  = halo_profile
        self.hod           = hod
        self.p_lin         = p_lin

        self.cosmo         = p_lin.cosmo

        self.rho_M0        = self.cosmo.omm * self.cosmo.H0**2 * 27751973.7

    def p_kz(self, k, z):
        return self.p_kz_1h(k, z) + self.p_kz_2h(k, z)

    def p_kz_1h(self, k, z):
        def integrand(logm):
            m = np.exp(logm)
            return m * self.mass_function.dndM_mz(m,z) * (self.hod.ngal(m, z) / self.hod.nbar(z)) * (m / self.rho_M0) * self.halo_profile.u_km(k,m,z)**2
        return scipy.integrate.quad( integrand, np.log(self.mass_function.Mmin), np.log(self.mass_function.Mmax) )[0]

    def p_kz_2h(self, k, z):
        def integrand_dg(logm):
            m = np.exp(logm)
            return m * self.mass_function.dndM_mz(m,z) * (self.hod.ngal(m,z) / self.hod.nbar(z)) * self.mass_function.b_mz(m, z) * self.halo_profile.u_km(k,m,z)
        # assume \int dM dN/dM b(M) M/\rho u(k,M) = 1
        
        return scipy.integrate.quad( integrand_dg, np.log(self.mass_function.Mmin), np.log(self.mass_function.Mmax) )[0] * self.p_lin.p_kz(k, z)
