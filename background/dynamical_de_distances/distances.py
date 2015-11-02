import numpy as np
np.seterr(all="ignore")
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import scipy.integrate
from distance_calc import *
cosmo = section_names.cosmological_parameters
dist = section_names.distances


def setup(options):
	section = option_section
	verbose = options.get_bool(section, "verbose", default=False)
	wmodel = options.get_int(section, "w_model", default=0)
	zmin = options.get_double(section, "zmin", default=0.0)
	zmax = options.get_double(section, "zmax", default=4.)
	dz = options.get_double(section, "dz", default=0.01)
	return zmin,zmax,dz,verbose,wmodel


def execute(block, config):
	zmin,zmax,dz,verbose,wmodel = config

	h0 = block[cosmo, 'h0']
	omega_m = block[cosmo, 'omega_m']
	omega_lambda = block[cosmo, 'omega_lambda']
	omega_k = block[cosmo, 'omega_k']
        if wmodel ==1:
	    de_params = (block[cosmo, 'w0'],block[cosmo, 'wa'],block[cosmo, 'ap'])
        elif wmodel ==2:
            de_params = (block[cosmo, 'w0'],block[cosmo, 'ode_e'])
        else:
	    de_params = (block[cosmo, 'w0'],block[cosmo, 'wa'])

        if verbose: print "om,ol,ok,h0,de_params", omega_m,omega_lambda, omega_k ,h0,de_params

        DM = DistanceCalc(omega_m,omega_k,omega_lambda,wmodel,de_params,h0)
        z_array = np.arange(zmin,zmax+dz,dz)
        mu = np.zeros_like(z_array)
        mu = map(DM.mu,z_array)
        da = map(DM.d_a,z_array)
        dl = map(DM.d_l,z_array)
             
        block[dist,'d_l'] = dl[::-1] # reverse order to match output from camb
        block[dist,'d_a'] = da[::-1]
	block[dist, 'mu'] = mu[::-1]
	block[dist, 'z'] = z_array[::-1]
	return 0

def cleanup(config):
    return 0
