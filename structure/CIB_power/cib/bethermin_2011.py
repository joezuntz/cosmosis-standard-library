# Number counts from Bethermin et. al. (2011) arxiv:1010.1150
# Note: requires scipy > 0.9.0, for IDLsave.

import os, subprocess
import numpy as np

import scipy.io
import scipy.interpolate

import quickspec as qs

basedir = os.environ.get('QUICKSPEC_DATA', os.path.dirname(qs.__file__) + "/data") + "/Bethermin_2011/"

class kern():
    def __init__(self, nu, smax, cosmo, counts):
        self.nu = nu
        self.smax = smax
        self.cosmo = cosmo
        self.counts = counts

    def w_lxz(self, l, x, z):
        return 1./(1.+z) * self.counts.jbar(self.nu, z, cosmo=self.cosmo, smax=self.smax)
    
class counts():
    def __init__(self, model="mean"):
        print "quickspec::cib::bethermin_2011::counts:: loading " + model + "mean"

        tfname = "dndsnudz_arr_" + model + "model_final.save"
        if not os.path.exists(basedir + tfname):
            if not os.path.exists(basedir): os.makedirs(basedir)
            qs.util.download("http://www.ias.u-psud.fr/irgalaxies/Model/save/" + tfname + ".gz", basedir + tfname + ".gz")
            subprocess.call(['gunzip', basedir + tfname + ".gz"])
        sav = scipy.io.idl.readsav( basedir + tfname )
        
        self.ls = sav['lambda']
        self.zs = sav['z']
        self.ss = sav['snu']
        self.dndsdz = sav['dndsnudz_arr'] # [ l, z, s ]

        self.ds = qs.util.deriv(self.ss)
        self.dz = qs.util.deriv(self.zs)

        # Note: Assign
        #    143 -> 2100 rather than 2097
        self.nu2ls = { 143.e9 : 2100, 150.e9 : 2000, 217.e9 : 1380, 220.e9 : 1360, 353.e9 : 850, 545.e9 : 550, 600.e9 : 500, 857.e9 : 350, 1200.e9 : 250 }

    def dNdS(self, nu, s, zmin=0., zmax=1100.):
        l = self.nu2ls[nu]
        il = np.where( self.ls == l )[0][0]

        si_fid = np.argmin( (self.ss - s)**2 )

        si_min = max(si_fid - 1, 0)
        si_max = min(si_min + 3, len(self.ss))
        si_min = max(si_max - 3, 0)

        ss = []
        rs = []
        for si in xrange(si_min, si_max):
            ss.append( self.ss[si] )

            zidx = np.where( (self.zs <= zmax) * (self.zs >= zmin) )[0]
        
            if len(zidx) == 0:
                rs.append( 0.0 )
            else:
                rs.append( np.sum( self.dz[zidx] *
                                   self.dndsdz[il,zidx,si] ) )
        
        return qs.interp.lagrange(s, np.array(ss), np.array(rs))

    def jbar(self, nu, z, smax=None, cosmo=None):
        # \bar{j}(nu, z) = (1+z) \int_0^{Scut} dS S [d^2N/dSdz]

        l = self.nu2ls[nu]
        il = np.where( self.ls == l )[0][0]
    
        iz_fid = np.argmin( (self.zs - z)**2 )

        iz_min = max(iz_fid - 1, 0)
        iz_max = min(iz_min + 3, len(self.zs))
        iz_min = max(iz_max - 3, 0)

        zs = []
        rs = []
        for iz in xrange(iz_min, iz_max):
            zs.append( self.zs[iz] )

            if smax != None:
                sidx = np.where( self.ss < smax )[0]
            else:
                sidx = np.arange(0, len(self.ss))
        
            if len(sidx) == 0:
                rs.append( 0.0 )
            else:
                rs.append( np.sum( self.ds[sidx] *
                                   self.ss[sidx] *
                                   self.dndsdz[il,iz,:][sidx] ) )
        
        return (1.+z) * qs.interp.lagrange(z, np.array(zs), np.array(rs)) * cosmo.H_z(z) / 3.e5

class jbar_pep():
    """ Bethermin 2011 jbar, precomputed from analytical formula by Penin
    NOTE: this jbar includes a color correction, while counts (above) does not. """
    def __init__(self):
        self.smaxs = {217e9 : 160.e-3, 353e9 : 325.e-3,  545e9 : 540.e-3, 857e9 : 710.e-3 }

        self.counts_spl = {}
        for frq in ['217', '353', '545', '857']:
            dat = np.loadtxt(os.path.dirname(__file__) + "/data/Bethermin_2011_jbar/j_z_" + frq + "GHz.dat")
            self.counts_spl[float(frq)*1.e9] = scipy.interpolate.UnivariateSpline( dat[:,0].flatten(), dat[:,1].flatten(), k=3, s=0 )

    def jbar(self, nu, z):
        return self.counts_spl[nu](z)
