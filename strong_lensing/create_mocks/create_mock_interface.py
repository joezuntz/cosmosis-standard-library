from cosmosis.datablock import names, option_section
from time_delay_likelihood import TimeDelayLikelihood, B1608, RXJ1131
from scipy.stats import lognorm
import numpy as np

def draw_from_lognorm(mu,sig,lbd):
    return lognorm.rvs(sig, loc=mu,size=1,scale=lbd)

def D_deltat(z_d,z_s, comovingDistance, omega_k, H0):
    c = 299792.4580  #km/s
    D_H = c/H0 #Mpc
    chi_s = comovingDistance(z_s)
    chi_d = comovingDistance(z_d)

    D_s = chi_s / (1+z_s)
    D_d = chi_d / (1+z_d)

    f_s = np.sqrt(1+omega_k*chi_d**2/D_H)
    f_d = np.sqrt(1+omega_k*chi_s**2/D_H)
    D_ds = (f_s*chi_s - f_d*chi_d) / (1+z_s)
    return (1+z_d) * D_d * D_s / D_ds


def setup(options):

	sig = options[option_section, "sigma"]
	lbd = options[option_section, "lambdaD"]
        fname = options[option_section, "filename"]
        data = np.loadtxt(fname)
        z_d = data[:,0]; z_s = data[:,1] 
	return sig,lbd, z_d,z_s


def execute(block, config):
	sig,lbd, z_d,z_s = config

	z_m = block[names.distances, "z"][::-1]
	d_m = block[names.distances, "d_m"][::-1]
	omega_k = block[names.cosmological_parameters, "omega_k"]
	H0 = block[names.cosmological_parameters, "hubble"]
	comovingDistance = lambda z: np.interp(z, z_m, d_m)

        Ddt_obs = np.zeros(len(z_d))
        for i,z in enumerate(z_d):
            Ddt_true=D_deltat(z,z_s[i], comovingDistance, omega_k, H0)
            d_obs=draw_from_lognorm(Ddt_true,sig,lbd)
            Ddt_obs[i]=d_obs
            #print d_obs, Ddt_true
        err=np.ones(len(Ddt_obs))*sig
        np.savetxt('cosmosis-standard-library/strong_lensing/time_delay_lenses/mock.txt',np.vstack((z_d,z_s,Ddt_obs,err)).T)
	return 0



