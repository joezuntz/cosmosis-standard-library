from builtins import object
import numpy as np
import scipy.integrate
c_km_per_s = 299792.458


class DistanceCalc(object):
    def __init__(self, om, ok, ol, wmodel, de_params, h0):
        self.om = om
        self.ol = ol
        self.ok = ok
        if self.ok == 0:
            self.is_flat = True
        self.wmodel = wmodel
        self.de_params = de_params
        self.h0 = h0
        self.d_h = c_km_per_s / (100.)  # Mpc/h

    def wfunc(self, a):
        if self.wmodel == 0:  # e.g. Linder 2003
            w0, wa = self.de_params
            return w0 + (1.0 - a) * wa
        elif self.wmodel == 1:  # e.g Huterer & Turner 2001
            w0, wa, ap = self.de_params
            return w0 + (ap - a) * wa
        elif self.wmodel == 2:  # Wetterich 2004
            w0, ode_e = self.de_params
            b = -3. * (w0 / (np.log((1. - ode_e) / ode_e) +
                             np.log((1. - self.om) / self.om)))
            return w0 / (1.0 - b * (np.log(a)))

    def w_integrand(self, a):
        return (1. + self.wfunc(a)) / a

    def w_int(self, z):
        a = 1. / (1. + z)
        return 3.0 * scipy.integrate.quad(self.w_integrand, a, 1.0, epsrel=1e-6, limit=50)[0]

    def e_z_inverse(self, z):
        return 1. / (np.sqrt(self.om * (1 + z)**3 + self.ok * (1 + z)**2 + self.ol * np.exp(self.w_int(z))))

    def d_c(self, z):
        return self.d_h * scipy.integrate.quad(self.e_z_inverse, 0.0, z, epsrel=1e-6, limit=50)[0]

    def d_m(self, z):
        '''returns los comoving distance in Mpc/h '''
        if self.is_flat:
            return self.d_c(z)
        else:
            chi = np.sqrt(np.abs(self.ok)) * self.d_c(z) / self.d_h
            if self.ok > 0:
                return self.d_h / np.sqrt(np.abs(self.ok)) * np.sinh(chi)
            elif self.ok < 0:
                return self.d_h / np.sqrt(np.abs(self.ok)) * np.sin(chi)

    def hubble(self, z):
        return 100. * self.h0 * 1. / self.e_z_inverse(z)

    def d_l(self, z):
        '''luminosity dist in Mpc '''
        return (1. + z) * self.d_m(z) / self.h0

    def d_a(self, z):
        '''Ang diameter distance in Mpc'''
        return self.d_m(z) / (1. + z) / self.h0

    def mu(self, z):
        return 5. * np.log10((1. + z) * self.d_m(z) / self.h0) + 25.
