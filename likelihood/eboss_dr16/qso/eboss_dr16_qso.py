import sys
import os

from numpy import log, pi, interp, where, loadtxt,dot, append, linalg

dirname = os.path.split(__file__)[0]

sys.path.append(dirname + '/../lrg')
from eboss_dr16_lrg import Eboss16LRGLikelihood

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]


class Eboss16QSOLikelihood(Eboss16LRGLikelihood):
    
    data_type = "QSO"
    like_name = "eboss16_qso"
    
    def build_data(self):
        
        self.mode = self.options.get_int("mode", default=0)
        
        if not self.mode:
            print('QSO data')
            print("BAO only: Dm(z)/rd and Dh(z)/rd")
            # Reading data file
            DATA_file = os.path.join(ROOT_dir, 'sdss_DR16_QSO_BAO_DMDH.txt')
        else:
            print('QSO data')
            print("BAO+FS: Dm(z)/rd, Dh(z)/rd and f(z)sigma8(z)")
            # Reading data file
            DATA_file = os.path.join(ROOT_dir, 'sdss_DR16_QSO_FSBAO_DMDHfs8.txt')
            
        DATA = loadtxt(DATA_file, usecols=(0,1))
        z_eff, data = DATA[:,0], DATA[:,1]
        
        return z_eff, data
    
    def build_covariance(self):
        
        if not self.mode:
            # Reading covariance matrix file
            COV_file = os.path.join(ROOT_dir, 'sdss_DR16_QSO_BAO_DMDH_covtot.txt')
        else:
            # Reading covariance matrix file
            COV_file = os.path.join(ROOT_dir, 'sdss_DR16_QSO_FSBAO_DMDHfs8_covtot.txt')
            
        cov = loadtxt(COV_file)
        self.inv_cov = linalg.inv(cov)
        return cov
    
setup, execute, cleanup = Eboss16QSOLikelihood.build_module()
