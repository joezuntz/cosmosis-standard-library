'''
Python version of Planck's plik-lite likelihood with the option to include
the low-ell temperature as two Gaussian bins

The official Planck likelihoods are availabe at https://pla.esac.esa.int/
The papers describing the Planck likelihoods are
Planck 2018: https://arxiv.org/abs/1907.12875
Planck 2015: https://arxiv.org/abs/1507.02704

The covariance matrix treatment is based on Zack Li's ACT likelihood code
available at: https://github.com/xzackli/actpols2_like_py

planck calibration is set to 1 by default but this can easily be modified
'''
import numpy as np
from scipy.io import FortranFile
import scipy.linalg

def main():
    TTTEEE2018=PlanckLitePy(year=2018, spectra='TTTEEE', use_low_ell_bins=False)
    TTTEEE2018.test()

    TTTEEE2018_lowTTbins=PlanckLitePy(year=2018, spectra='TTTEEE', use_low_ell_bins=True)
    TTTEEE2018_lowTTbins.test()

    TT2018=PlanckLitePy(year=2018, spectra='TT', use_low_ell_bins=False)
    TT2018.test()

    TT2018_lowTTbins=PlanckLitePy(year=2018, spectra='TT', use_low_ell_bins=True)
    TT2018_lowTTbins.test()


class PlanckLitePy:
    def __init__(self, data_directory='data', year=2018, spectra='TT', use_low_ell_bins=False):
        '''
        data_directory = path from where you are running this to the folder
          containing the planck2015/8_low_ell and planck2015/8_plik_lite data
        year = 2015 or 2018
        spectra = TT for just temperature or TTTEEE for temperature (TT),
          E mode (EE) and cross (TE) spectra
        use_low_ell_bins = True to use 2 low ell bins for the TT 2<=ell<30 data
          or False to only use ell>=30
        '''
        self.year=year
        self.spectra=spectra
        self.use_low_ell_bins=use_low_ell_bins #False matches Plik_lite - just l>=30

        if self.use_low_ell_bins:
            self.nbintt_low_ell=2
            self.plmin_TT=2
        else:
            self.nbintt_low_ell=0
            self.plmin_TT=30
        self.plmin=30
        self.plmax=2508
        self.calPlanck=1

        if year==2015:
            self.data_dir=data_directory+'/planck2015_plik_lite/'
            version=18
        elif year==2018:
            self.data_dir=data_directory+'/planck2018_plik_lite/'
            version=22
        else:
            print('Year must be 2015 or 2018')
            return 1

        if spectra=='TT':
            self.use_tt=True
            self.use_ee=False
            self.use_te=False
        elif spectra=='TTTEEE':
            self.use_tt=True
            self.use_ee=True
            self.use_te=True
        else:
            print('Spectra must be TT or TTTEEE')
            return 1

        self.nbintt_hi = 215 #30-2508   #used when getting covariance matrix
        self.nbinte = 199 #30-1996
        self.nbinee = 199 #30-1996
        self.nbin_hi=self.nbintt_hi+self.nbinte+self.nbinee

        self.nbintt=self.nbintt_hi+self.nbintt_low_ell #mostly want this if using low ell
        self.nbin_tot=self.nbintt+self.nbinte+self.nbinee

        self.like_file = self.data_dir+'cl_cmb_plik_v'+str(version)+'.dat'
        self.cov_file  = self.data_dir+'c_matrix_plik_v'+str(version)+'.dat'
        self.blmin_file = self.data_dir+'blmin.dat'
        self.blmax_file = self.data_dir+'blmax.dat'
        self.binw_file = self.data_dir+'bweight.dat'

        # read in binned ell value, C(l) TT, TE and EE and errors
        # use_tt etc to select relevant parts
        self.bval, self.X_data, self.X_sig=np.genfromtxt(self.like_file, unpack=True)
        self.blmin=np.loadtxt(self.blmin_file).astype(int)
        self.blmax=np.loadtxt(self.blmax_file).astype(int)
        self.bin_w=np.loadtxt(self.binw_file)

        if self.use_low_ell_bins:
            self.data_dir_low_ell=data_directory+'/planck'+str(year)+'_low_ell/'
            self.bval_low_ell, self.X_data_low_ell, self.X_sig_low_ell=np.genfromtxt(self.data_dir_low_ell+'CTT_bin_low_ell_'+str(year)+'.dat', unpack=True)
            self.blmin_low_ell=np.loadtxt(self.data_dir_low_ell+'blmin_low_ell.dat').astype(int)
            self.blmax_low_ell=np.loadtxt(self.data_dir_low_ell+'blmax_low_ell.dat').astype(int)
            self.bin_w_low_ell=np.loadtxt(self.data_dir_low_ell+'bweight_low_ell.dat')

            self.bval=np.concatenate((self.bval_low_ell, self.bval))
            self.X_data=np.concatenate((self.X_data_low_ell, self.X_data))
            self.X_sig=np.concatenate((self.X_sig_low_ell, self.X_sig))

            self.blmin_TT=np.concatenate((self.blmin_low_ell, self.blmin+len(self.bin_w_low_ell)))
            self.blmax_TT=np.concatenate((self.blmax_low_ell, self.blmax+len(self.bin_w_low_ell)))
            self.bin_w_TT=np.concatenate((self.bin_w_low_ell, self.bin_w))

        else:
            self.blmin_TT=self.blmin
            self.blmax_TT=self.blmax
            self.bin_w_TT=self.bin_w

        self.mu = self._cut_vector(self.X_data)
        self.fisher=self.get_inverse_covmat()

    def get_inverse_covmat(self):
        #read full covmat
        f = FortranFile(self.cov_file, 'r')
        covmat = f.read_reals(dtype=float).reshape((self.nbin_hi,self.nbin_hi))
        for i in range(self.nbin_hi):
            for j in range(i,self.nbin_hi):
                covmat[i,j] = covmat[j,i]

        #select relevant covmat
        if self.use_tt and not(self.use_ee) and not(self.use_te):
            #just tt
            bin_no=self.nbintt_hi
            start=0
            end=start+bin_no
            cov=covmat[start:end, start:end]
        elif not(self.use_tt) and not(self.use_ee) and self.use_te:
            #just te
            bin_no=self.nbinte
            start=self.nbintt_hi
            end=start+bin_no
            cov=covmat[start:end, start:end]
        elif not(self.use_tt) and self.use_ee and not(self.use_te):
            #just ee
            bin_no=self.nbinee
            start=self.nbintt_hi+self.nbinte
            end=start+bin_no
            cov=covmat[start:end, start:end]
        elif self.use_tt and self.use_ee and self.use_te:
            #use all
            bin_no=self.nbin_hi
            cov=covmat
        else:
            raise NotImplementedError("Specified combination of tt, te, ee is not implemented")

        #invert high ell covariance matrix (cholesky decomposition should be faster)
        fisher=scipy.linalg.cho_solve(scipy.linalg.cho_factor(cov), np.identity(bin_no))
        fisher=fisher.transpose()


        if self.use_low_ell_bins:
            bin_no += self.nbintt_low_ell
            inv_covmat_with_lo=np.zeros(shape=(bin_no, bin_no))
            inv_covmat_with_lo[0:2, 0:2]=np.diag(1./self.X_sig_low_ell**2)
            inv_covmat_with_lo[2:,2:]= fisher
            fisher=inv_covmat_with_lo
            fullcov=np.zeros(shape=(bin_no, bin_no))
            fullcov[0:2,0:2] = np.diag(self.X_sig_low_ell**2)
            fullcov[2:,2:] = cov
            cov = fullcov

        self.cov = cov

        return fisher

    def make_mean_vector(self, Dltt, Dlte, Dlee, ellmin=2):
        #convert model Dl's to Cls then bin them
        ls=np.arange(len(Dltt))+ellmin
        fac=ls*(ls+1)/(2*np.pi)
        Cltt=Dltt/fac
        Clte=Dlte/fac
        Clee=Dlee/fac


        # Fortran to python slicing: a:b becomes a-1:b
        # need to subtract 1 to use 0 indexing for cl,
        # then add one for weights because fortran includes top value
        Cltt_bin=np.zeros(self.nbintt)
        for i in range(self.nbintt):
            Cltt_bin[i]=np.sum(Cltt[self.blmin_TT[i]+self.plmin_TT-ellmin:self.blmax_TT[i]+self.plmin_TT+1-ellmin]*self.bin_w_TT[self.blmin_TT[i]:self.blmax_TT[i]+1])

        # bin widths and weights are the same for TT, TE and EE
        Clte_bin=np.zeros(self.nbinte)
        for i in range(self.nbinte):
            Clte_bin[i]=np.sum(Clte[self.blmin[i]+self.plmin-ellmin:self.blmax[i]+self.plmin+1-ellmin]*self.bin_w[self.blmin[i]:self.blmax[i]+1])

        # bin widths and weights are the same for TT, TE and EE
        Clee_bin=np.zeros(self.nbinee)
        for i in range(self.nbinee):
            Clee_bin[i]=np.sum(Clee[self.blmin[i]+self.plmin-ellmin:self.blmax[i]+self.plmin+1-ellmin]*self.bin_w[self.blmin[i]:self.blmax[i]+1])

        X_model=np.zeros(self.nbin_tot)
        X_model[:self.nbintt]=Cltt_bin/self.calPlanck**2
        X_model[self.nbintt:self.nbintt+self.nbinte]=Clte_bin/self.calPlanck**2
        X_model[self.nbintt+self.nbinte:]=Clee_bin/self.calPlanck**2

        return X_model

    def _cut_vector(self, Y):
        #choose relevant bits based on whether using TT, TE, EE
        if self.use_tt and not(self.use_ee) and not(self.use_te):
            #just tt
            bin_no=self.nbintt
            start=0
            end=start+bin_no
            v=Y[start:end]
        elif not(self.use_tt) and not(self.use_ee) and self.use_te:
            #just te
            bin_no=self.nbinte
            start=self.nbintt
            end=start+bin_no
            v=Y[start:end]
        elif not(self.use_tt) and self.use_ee and not(self.use_te):
            #just ee
            bin_no=self.nbinee
            start=self.nbintt+self.nbinte
            end=start+bin_no
            v=Y[start:end]
        elif self.use_tt and self.use_ee and self.use_te:
            #use all
            bin_no=self.nbin_tot
            v=Y
        else:
            raise NotImplementedError("Specified combination of tt, te, ee is not implemented")
        return v

    def loglike(self, Dltt, Dlte, Dlee, ellmin=2):

        X_model = self.make_mean_vector(Dltt, Dlte, Dlee, ellmin=ellmin)

        Y=self.X_data - X_model

        diff_vec = self._cut_vector(Y)

        return -0.5*diff_vec.dot(self.fisher.dot(diff_vec))


    def test(self):
        ls, Dltt, Dlte, Dlee = np.genfromtxt('data/Dl_planck2015fit.dat', unpack=True)
        ellmin=int(ls[0])
        loglikelihood=self.loglike(Dltt, Dlte, Dlee, ellmin)

        if self.year==2018 and self.spectra=='TTTEEE' and not self.use_low_ell_bins:
            print('Log likelihood for 2018 high-l TT, TE and EE:')
            expected = -291.33481235418026
            # Plik-lite within cobaya gives  -291.33481235418003
        elif self.year==2018 and self.spectra=='TTTEEE' and self.use_low_ell_bins:
            print('Log likelihood for 2018 high-l TT, TE and EE + low-l TT bins:')
            expected = -293.95586501795134
        elif self.year==2018 and self.spectra=='TT' and not self.use_low_ell_bins:
            print('Log likelihood for 2018 high-l TT:')
            expected = -101.58123068722583
            #Plik-lite within cobaya gives -101.58123068722568
        elif self.year==2018 and self.spectra=='TT' and self.use_low_ell_bins:
            print('Log likelihood for 2018 high-l TT + low-l TT bins:')
            expected = -104.20228335099686

        elif self.year==2015 and self.spectra=='TTTEEE' and not self.use_low_ell_bins:
            print('NB: Don\'t use 2015 polarization!')
            print('Log likelihood for 2015 high-l TT, TE and EE:')
            expected = -280.9388125627618
            # Plik-lite within cobaya gives  -291.33481235418003
        elif self.year==2015 and self.spectra=='TTTEEE' and self.use_low_ell_bins:
            print('NB: Don\'t use 2015 polarization!')
            print('Log likelihood for 2015 high-l TT, TE and EE + low-l TT bins:')
            expected = -283.1905700256343
        elif self.year==2015 and self.spectra=='TT' and not self.use_low_ell_bins:
            print('Log likelihood for 2015 high-l TT:')
            expected = -102.34403873289027
            #Plik-lite within cobaya gives -101.58123068722568
        elif self.year==2015 and self.spectra=='TT' and self.use_low_ell_bins:
            print('Log likelihood for 2015 high-l TT + low-l TT bins:')
            expected = -104.59579619576277
        else:
            expected=None

        print('Planck-lite-py:',loglikelihood)
        if(expected):
            print('expected:', expected)
            print('difference:', loglikelihood-expected, '\n')



if __name__=='__main__':
    main()
