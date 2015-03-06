import numpy as np

import scipy.interpolate

import quickspec as qs
import mps

class mps_lin_camb(mps.mps_lin):
    def __init__(self, cosmo, sips, kmax=200., logk_spacing=0.02, nonlinear=False,
                 zvec = np.array([2048., 1024., 512., 256., 128., 64.,
                                  32., 16., 12., 8., 6., 4., 3., 2.,
                                  1.5, 1., 0.5, 0.])):
        """ wrapper class for the linear matter power spectrum computed by CAMB. """

        try:
            import pycamb
        except:
            print "pycamb (https://github.com/joezuntz/pycamb) could not be loaded. is it installed?"
            raise

        self.cosmo     = cosmo

        par = {}
        par['H0']             = cosmo.H0
        par['omegab']         = cosmo.omb
        par['omegac']         = cosmo.omc
        par['omegav']         = cosmo.oml

        par['scalar_amp']     = sips.amp
        par['scalar_index']   = sips.n_s
        par['scalar_running'] = sips.n_r
        par['scalar_pivot']   = sips.k_pivot
        par['NonLinear']      = nonlinear

        par['maxk']           = kmax
        par['logk_spacing']   = logk_spacing
        self.par              = par

        p_kz = pycamb.matter_power( zvec, **par )

        self.arr_z = zvec[::-1]
        self.arr_k = p_kz[0][:,0] * (self.cosmo.h)
        self.mat_p = p_kz[1][:,::-1] / (self.cosmo.h)**3
        for iz, z in enumerate(self.arr_z):
            self.mat_p[:,iz] *= (1.+z)**2

        self.spl_p = scipy.interpolate.RectBivariateSpline(np.log(self.arr_k), self.arr_z, self.mat_p, kx=3, ky=3, s=0)

        self.kmin = self.arr_k[+0]*0.999
        self.kmax = self.arr_k[-1]*1.001

        self.zmin = np.min(zvec)
        self.zmax = np.max(zvec)

    def p_kz(self, k, z):
        """ returns the amplitude of the matter power spectrum at wavenumber k (in Mpc^{-1}) and conformal distance x (in Mpc). """

        assert( np.all( k >= self.kmin ) )
        assert( np.all( k <= self.kmax ) )
        assert( np.all( z >= self.zmin ) )
        assert( np.all( z <= self.zmax ) )

        k, z, s = qs.util.pair(k, z)
        ret  = self.spl_p.ev( np.log(k), z)
        ret /= (1.+z)**2
        return ret.reshape(s)
