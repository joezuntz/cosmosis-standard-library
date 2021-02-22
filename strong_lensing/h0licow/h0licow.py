# This likelihood incorporates the time-delay distances Dt
# and deflector distances Dd of 6 strongly lensed quasars
# from H0LiCOW into MontePython.
# Implementation by Sherry Suyu and Stefan Taubenberger.
# Version 1.0
# JZ: adapted for cosmosis

import os
import numpy as np
from math import log, sqrt, pi
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from cosmosis.datablock import names, SectionOptions
import urllib


DEFAULT_DATA_DIR = os.path.join(os.path.split(__file__)[0], "timedelay_6lenses")


# url for downloading data if not present
url_base = "https://raw.githubusercontent.com/shsuyu/H0LiCOW-public/master/MontePython_cosmo_sampling/data/timedelay_6lenses"


class TimeDelaySixLens:


    # initialization routine

    def __init__(self, options):
        super(TimeDelaySixLens, self).__init__()

        # specify which distances to include for inference (1=include, 0=exclude)
        self.data_dt_b1608 = options.get_int("data_dt_b1608", default=0)
        self.data_dd_dt_b1608 = options.get_int("data_dd_dt_b1608", default=1)
        self.data_dt_j1206 = options.get_int("data_dt_j1206", default=0)
        self.data_dd_dt_j1206 = options.get_int("data_dd_dt_j1206", default=1)
        self.data_dt_wfi2033 = options.get_int("data_dt_wfi2033", default=1)
        self.data_dt_he0435 = options.get_int("data_dt_he0435", default=1)
        self.data_dt_pg1115 = options.get_int("data_dt_pg1115", default=0)
        self.data_dd_dt_pg1115 = options.get_int("data_dd_dt_pg1115", default=1)
        self.data_dt_rxj1131 = options.get_int("data_dt_rxj1131", default=0)
        self.data_dd_dt_rxj1131 = options.get_int("data_dd_dt_rxj1131", default=1)
        self.data_directory = options.get_string("data_dir", default=DEFAULT_DATA_DIR)

        print("Loading H0LiCOW data from", self.data_directory)
        

        if (self.data_dd_dt_b1608 == 1 and self.data_dt_b1608 == 1):
            print('Warning: configuration conflict. data_dt_b1608 cannot be included in addition to data_dd_dt_b1608. Only data_dd_dt_b1608 used.')
        if (self.data_dd_dt_j1206 == 1 and self.data_dt_j1206 == 1):
            print('Warning: configuration conflict. data_dt_j1206 cannot be included in addition to data_dd_dt_j1206. Only data_dd_dt_j1206 used.')
        if (self.data_dd_dt_pg1115 == 1 and self.data_dt_pg1115 == 1):
            print('Warning: configuration conflict. data_dt_pg1115 cannot be included in addition to data_dd_dt_pg1115. Only data_dd_dt_pg1115 used.')
        if (self.data_dd_dt_rxj1131 == 1 and self.data_dt_rxj1131 == 1):
            print('Warning: configuration conflict. data_dt_rxj1131 cannot be included in addition to data_dd_dt_rxj1131. Only data_dd_dt_rxj1131 used.')
        
        missing = []
        for filename in [
            'J1206_final.csv',
            'wfi2033_dt_bic.dat',
            'HE0435_Ddt_AO+HST.dat',
            'PG1115_AO+HST_Dd_Ddt.dat',
            'RXJ1131_AO+HST_Dd_Ddt.dat',
            'B1608_Dd_Ddt_params.dat']:
            if not os.path.exists(os.path.join(self.data_directory,filename)):
                missing.append(filename)

        if missing:
            msg = "\n".join([
                "curl -o {}/{}  {}/{}".format(self.data_directory, f, url_base, urllib.parse.quote(f)) for f in missing])
            raise IOError("\n\nPlease download the H0LiCOW data files like this to use it:\n" + msg + "\n\n")
        
        # J1206
        # =====
        # read in the J1206 Dt and Dd distribution
        j1206chain = np.loadtxt(os.path.join(self.data_directory,'J1206_final.csv'), skiprows=1, delimiter=',')
        j1206dist = j1206chain.transpose() # j1206dist[0] is the Dt distribution, and j1206dist[1] is the Dd distribution
        self.j1206_dt_fit_kde = gaussian_kde(j1206dist[0])
        j1206_dd_dt_values = np.vstack([j1206dist[0],j1206dist[1]])
        self.j1206_dd_dt_fit_kde = gaussian_kde(j1206_dd_dt_values)

        # WFI2033
        # =======
        # read in the WFI2033 Dt distribution
        wfi2033chain = np.loadtxt(os.path.join(self.data_directory,'wfi2033_dt_bic.dat'), skiprows=1, delimiter=',')
        wfi2033dist = wfi2033chain.transpose() #wfi2033dist[0] is the Dt distribution, and wfi2033dist[1] is the corresponding weight
        self.wfi2033_dt_fit_kde = gaussian_kde(wfi2033dist[0], weights = wfi2033dist[1])
        
        # HE0435
        # ======
        # read in the HE0435 Dt distribution
        he0435chain = np.loadtxt(os.path.join(self.data_directory,'HE0435_Ddt_AO+HST.dat'))
        he0435dist = he0435chain
        self.he0435_dt_fit_kde = gaussian_kde(he0435dist)
        
        # PG1115
        # ======
        # read in the PG1115 Dt and Dd distribution
        pg1115chain = np.loadtxt(os.path.join(self.data_directory,'PG1115_AO+HST_Dd_Ddt.dat'), skiprows=1)
        pg1115dist = pg1115chain.transpose() # pg1115dist[0] is the Dd distribution, and pg1115dist[1] is the Dt distribution (note, opposite to J1206)
        self.pg1115_dt_fit_kde = gaussian_kde(pg1115dist[1]) 
        pg1115_dd_dt_values = np.vstack([pg1115dist[1],pg1115dist[0]])  #[Dt,Dd]
        self.pg1115_dd_dt_fit_kde = gaussian_kde(pg1115_dd_dt_values)

        # RXJ1131
        # =======
        # read in the RXJ1131 Dt and Dd distribution
        rxj1131chain = np.loadtxt(os.path.join(self.data_directory,'RXJ1131_AO+HST_Dd_Ddt.dat'), skiprows=1)
        rxj1131dist = rxj1131chain.transpose() # rxj1131dist[0] is the Dd distribution, and rxj1131dist[1] is the Dt distribution (note, opposite to J1206)
        self.rxj1131_dt_fit_kde = gaussian_kde(rxj1131dist[1]) 
        rxj1131_dd_dt_values = np.vstack([rxj1131dist[1],rxj1131dist[0]])  #[Dt,Dd]
        self.rxj1131_dd_dt_fit_kde = gaussian_kde(rxj1131_dd_dt_values)

        # B1608
        # =====
        # use analytic skewed log normal fit to Dd and Dt, so no need to read chain here and do kde
        b1608params = np.loadtxt(os.path.join(self.data_directory,'B1608_Dd_Ddt_params.dat'), skiprows=9)
        self.lambda_d = b1608params[2]
        self.mu_d = b1608params[3]
        self.sigma_d = b1608params[4]
        self.lambda_sft = b1608params[5]
        self.mu_sft = b1608params[6]
        self.sigma_sft = b1608params[7]



    # compute Dds given Dd and Ds
    
    def calculate_Dds (self, ok, K, zd, zs, Dd, Ds):

        # to compute Dds from Dd and Ds, first need to figure out whether the universe is flat
        
        if (np.fabs(ok)<1.e-6): #flat Universe, so can use simple formula            
            Dds = ((1. + zs) * Ds - (1 + zd) * Dd) / (1. + zs)

        elif (K > 0): 
            chis = np.arcsin ( Ds * (1. + zs) * np.sqrt(K) ) / np.sqrt(K)
            chid = np.arcsin ( Dd * (1. + zd) * np.sqrt(K) ) / np.sqrt(K)
            chids = chis - chid
            Dds = (1./(1+zs)) * (1./np.sqrt(K)) * np.sin (np.sqrt(K)*chids)
            
        else: #K<0
            chis = np.arcsinh ( Ds * (1. + zs) * np.sqrt(-K) ) / np.sqrt(-K)
            chid = np.arcsinh ( Dd * (1. + zd) * np.sqrt(-K) ) / np.sqrt(-K)
            chids = chis - chid
            Dds = (1./(1+zs)) * (1./np.sqrt(-K)) * np.sinh (np.sqrt(-K)*chids)

        return Dds

        
    # compute likelihood

    def likelihood(self, block):

        lkl = 0.

        # current cosmo parameter values
        
        h = block['cosmological_parameters', 'h0']
        ok = block['cosmological_parameters', 'omega_k']

        z_bg = block['distances', 'z']
        D_A = block['distances', 'D_A']
        angular_distance = interp1d(z_bg, D_A, kind='linear')


        c = 299792.458 #in km/s
        K = -1. * ok * h**2 * 10000. / (c**2)  # K is now in Mpc^-2
        
        if (self.data_dt_b1608 == 1 and self.data_dd_dt_b1608 == 0):
            # the following is to add the Dt skewed log normal likelihood of Suyu et al. (2010):
            
            zd1608 = 0.6304
            zs1608 = 1.394
            
            Dd = angular_distance(zd1608)
            Ds = angular_distance(zs1608)

            # to compute Dds from Dd and Ds
            Dds = self.calculate_Dds (ok, K, zd1608, zs1608, Dd, Ds)
            
            Dt = (1 + zd1608) * Dd * Ds / Dds
 
            if (Dt > self.lambda_d):
                lkl1608 = - (log(Dt - self.lambda_d) - self.mu_d) ** 2 / 2. / self.sigma_d ** 2 - log(sqrt(2. * pi) * (Dt - self.lambda_d) * self.sigma_d)
                
                lkl = lkl + lkl1608
                
            else:
                return -np.inf



        if (self.data_dd_dt_b1608 == 1):
            # the following is to add the Dt skewed log normal likelihood of Suyu et al. (2010) and the Dd skewed log normal likelihood of Jee et al. (2019):
            
            zd1608 = 0.6304
            zs1608 = 1.394

            Dd = angular_distance(zd1608)
            Ds = angular_distance(zs1608)

            # to compute Dds from Dd and Ds
            Dds = self.calculate_Dds (ok, K, zd1608, zs1608, Dd, Ds)
                
            Dt = (1 + zd1608) * Dd * Ds / Dds
            
            if (Dt > self.lambda_d and Dd > self.lambda_sft):
                lkl1608 = - (log(Dt - self.lambda_d) - self.mu_d) ** 2 / 2. / self.sigma_d ** 2 - log(sqrt(2. * pi) * (Dt - self.lambda_d) * self.sigma_d)
                lkl1608bothD = lkl1608 - (log(Dd - self.lambda_sft) - log(self.mu_sft)) ** 2 / 2. / self.sigma_sft ** 2 - log(sqrt(2. * pi) * (Dd - self.lambda_sft) * self.sigma_sft)
                
                lkl = lkl + lkl1608bothD
                
            else:
                return -np.inf


                    
        if (self.data_dt_j1206 == 1 and self.data_dd_dt_j1206 == 0):
            # the following is to add the Dt KDE likelihood of Birrer et al. (2019):
            
            zd1206 = 0.745
            zs1206 = 1.789
            
            Dd = angular_distance(zd1206)
            Ds = angular_distance(zs1206)
            
            Dds = self.calculate_Dds (ok, K, zd1206, zs1206, Dd, Ds)
            
            Dt = (1 + zd1206) * Dd * Ds / Dds
            
            lkl1206 = np.log(self.j1206_dt_fit_kde(Dt))
            
            lkl = lkl + lkl1206
                


        if (self.data_dd_dt_j1206 == 1):
            # the following is to add the (Dt,Dd) KDE likelihood of Birrer et al. (2019):

            zd1206 = 0.745
            zs1206 = 1.789
        
            Dd = angular_distance(zd1206)
            Ds = angular_distance(zs1206)

            Dds = self.calculate_Dds (ok, K, zd1206, zs1206, Dd, Ds)
                               
            Dt = (1 + zd1206) * Dd * Ds / Dds
            
            lkl1206bothD = np.log(self.j1206_dd_dt_fit_kde([Dt,Dd]))

            lkl = lkl + lkl1206bothD



        if (self.data_dt_wfi2033 == 1):
            # the following is to add the Dt KDE likelihood of Rusu et al. (2019):

            zd2033 = 0.6575
            zs2033 = 1.662

            Dd =  angular_distance(zd2033)
            Ds = angular_distance(zs2033)

            Dds = self.calculate_Dds (ok, K, zd2033, zs2033, Dd, Ds)

            Dt = (1 + zd2033) * Dd * Ds / Dds

            lkl2033 = np.log(self.wfi2033_dt_fit_kde(Dt))
            
            lkl = lkl + lkl2033

 
 
        if (self.data_dt_he0435 == 1):
            # the following is to add the Dt KDE likelihood of Chen et al. (2019):

            zd0435 = 0.4546
            zs0435 = 1.693
        
            Dd = angular_distance(zd0435)
            Ds = angular_distance(zs0435)

            Dds = self.calculate_Dds (ok, K, zd0435, zs0435, Dd, Ds)
                               
            Dt = (1 + zd0435) * Dd * Ds / Dds
            
            lkl0435 = np.log(self.he0435_dt_fit_kde(Dt))

            lkl = lkl + lkl0435
           
           

        if (self.data_dt_pg1115 == 1 and self.data_dd_dt_pg1115 == 0):
            # the following is to add the Dt KDE likelihood of Chen et al. (2019):

            zd1115 = 0.311
            zs1115 = 1.722
            
            Dd = angular_distance(zd1115)
            Ds = angular_distance(zs1115)
            
            Dds = self.calculate_Dds (ok, K, zd1115, zs1115, Dd, Ds)
                       
            Dt = (1 + zd1115) * Dd * Ds / Dds
            
            lkl1115 = np.log(self.pg1115_dt_fit_kde(Dt))
            
            lkl = lkl + lkl1115

            
            
        if (self.data_dd_dt_pg1115 == 1):
            # the following is to add the (Dt,Dd) KDE likelihood of Chen et al. (2019):

            zd1115 = 0.311
            zs1115 = 1.722
        
            Dd = angular_distance(zd1115)
            Ds = angular_distance(zs1115)

            Dds = self.calculate_Dds (ok, K, zd1115, zs1115, Dd, Ds)
                               
            Dt = (1 + zd1115) * Dd * Ds / Dds
            
            lkl1115bothD = np.log(self.pg1115_dd_dt_fit_kde([Dt,Dd]))

            lkl = lkl + lkl1115bothD


                    
        if (self.data_dt_rxj1131 == 1 and self.data_dd_dt_rxj1131 == 0):
            # the following is to add the Dt KDE likelihood of Chen et al. (2019):
            
            zd1131 = 0.295
            zs1131 = 0.654
            
            Dd = angular_distance(zd1131)
            Ds = angular_distance(zs1131)
            
            Dds = self.calculate_Dds (ok, K, zd1131, zs1131, Dd, Ds)
                           
            Dt = (1 + zd1131) * Dd * Ds / Dds
            
            lkl1131 = np.log(self.rxj1131_dt_fit_kde(Dt))
            
            lkl = lkl + lkl1131
                
                

        if (self.data_dd_dt_rxj1131 == 1):
            # the following is to add the (Dt,Dd) KDE likelihood of Chen et al. (2019):

            zd1131 = 0.295
            zs1131 = 0.654
        
            Dd = angular_distance(zd1131)
            Ds = angular_distance(zs1131)

            Dds = self.calculate_Dds (ok, K, zd1131, zs1131, Dd, Ds)
                               
            Dt = (1 + zd1131) * Dd * Ds / Dds
            
            lkl1131bothD = np.log(self.rxj1131_dd_dt_fit_kde([Dt,Dd]))

            lkl = lkl + lkl1131bothD


                        
        return float(lkl)



def setup(options):
    options = SectionOptions(options)
    return TimeDelaySixLens(options)

def execute(block, config):
    like = config.likelihood(block)
    block[names.likelihoods, "h0licow_like"] = like
    return 0
