# CIB number counts from Lagache, Dole, Puget (2004), astro-ph/0209115
# Note: requires scipy > 0.9.0, for IDLsave.

import os
import numpy as np

import scipy.io

import quickspec as qs

basedir = os.environ.get('QUICKSPEC_DATA', os.path.dirname(qs.__file__) + "/data") + "/LDP_2004/"

class counts(object):
    """ CIB number counts from Lagache, Dole, Puget (2004), astro-ph/0209115 """
    def __init__(self, cold=False):
        self.nu2ls = { 143.e9 : 2097, 217.e9 : 1380, 353.e9 : 850, 545.e9 : 550, 857.e9 : 350, 1200.e9 : 250 }
            
        self.cold   = cold
        self.loaded = None

    def load(self, nu):
        if self.loaded == nu:
            return

        print "quickspec::cib::ldp_2004::counts:: loading " + str(nu/1e9) + "GHz"
        assert( nu in self.nu2ls.keys() )

        tfname = "create_counts_" + ('%04d' % self.nu2ls[nu]) + "_Omega_lambda" + [".save", ".cold.save"][self.cold]

        if not os.path.exists(basedir + tfname):
            if not os.path.exists(basedir): os.makedirs(basedir)
            qs.util.download("http://www.ias.u-psud.fr/irgalaxies/Model/save/" + tfname, basedir + tfname)
        
        sav = scipy.io.idl.readsav(basedir + tfname)

        self.zs   = sav['z']          # redshift
        self.ls   = sav['lum_array']  # luminosity L
        self.slz  = sav['slz']        # flux in Jy [z, L]

        self.dz   = qs.util.deriv( self.zs )
        self.dlnl = qs.util.deriv( np.log(self.ls) )

        self.dndlnldz = sav['dndlnldz']
        self.dslz = np.zeros( np.shape(self.slz) )
        for iz, z in enumerate(self.zs):
            self.dslz[iz, :] = qs.util.deriv( self.slz[iz, :] )

        self.loaded = nu

    def dNdS(self, nu, s, zmin=0., zmax=1100.):
        self.load(nu)

        ret = 0.0
        for iz in np.argwhere( (self.zs <= zmax) * (self.zs >= zmin) ):
            iz = iz[0]

            si_fid = np.argmin( (self.slz[iz,:] - s)**2 )

            si_min = max(si_fid - 1, 0)
            si_max = min(si_min + 3, len(self.ls))
            si_min = max(si_max - 3, 0)

            if (self.slz[iz,si_min] > s) or (self.slz[iz,si_max-1] < s):
                continue

            ret += self.dz[iz] * qs.interp.lagrange(s, self.slz[iz,si_min:si_max], self.dndlnldz[iz,si_min:si_max] * self.dlnl[si_min:si_max] / self.dslz[iz, si_min:si_max])

        return ret

    def jbar(self, nu, z, cosmo, smax=None):
        # \bar{j}(nu, z) = (1+z) \int_0^{Smax} dS S [d^2N/dSdz]
        #                = (1+z) \int dln(L) S(L,z) Heaviside(S-Smax) d^2N/dln(L)dz * H(z)
        #
        self.load(nu)
        
        iz_fid = np.argmin( (self.zs - z)**2 )

        iz_min = max(iz_fid - 1, 0)
        iz_max = min(iz_min + 3, len(self.zs))
        iz_min = max(iz_max - 3, 0)

        zs = []
        rs = []
        for iz in xrange(iz_min, iz_max):
            zs.append( self.zs[iz] )

            if smax != None:
                sidx = np.where( self.slz[iz,:] < smax )[0]
            else:
                sidx = np.arange(0, len(self.slz[iz,:]))
        
            if len(sidx) == 0:
                rs.append( 0.0 )
            else:
                rs.append( np.sum( self.dlnl[sidx] *
                                   self.slz[iz,sidx] *
                                   self.dndlnldz[iz,sidx] ) )

        return (1.+z) * qs.interp.lagrange(z, np.array(zs), np.array(rs)) * cosmo.H_z(z) / 3.e5
