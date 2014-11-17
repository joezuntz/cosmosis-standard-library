import numpy as np
import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import scipy.integrate
from math import log10, floor
from scipy.interpolate import UnivariateSpline, RectBivariateSpline,SmoothBivariateSpline

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
evs = section_names.evs
mf  = section_names.mass_function
dist = section_names.distances


def setup(options):
	section = option_section
	feedback = options.get_int(section, "feedback", default=0)
	redshift = options.get_double(section, "redshift", default=1.)
	output_pdf = options.get_int(section, "output_pdf", default=0)
	frac = options.get_double(section, "frac", default=1.0) #
	Mmin = options.get_double(section, "Mmin", default=1.e14) #
	Mmax= options.get_double(section, "Mmax", default=5.e15) #
	dm = options.get_double(section, "dm", default=5.e12) #
	minput = np.arange(Mmin,Mmax,dm)
	return (feedback, redshift,frac,minput,output_pdf)


def massfunction(m,zz,rbs):
	z=round(zz,2)
	return (rbs.ev(m,z)[0]) # This spline can be produce Nan or negative values for fm below if the mf has been saved on a coarse grid. soln is to increase nsteps in massfunction module

def dndmint(logm,zz,rbs):
	m = np.exp(logm)
	integrand = massfunction(m,zz,rbs)
	return integrand


def dvdm_zint( zz,m,omega_matter,h0,interp_da,rbs):
	return dVcdz(zz,omega_matter,h0,interp_da) *(1./m)*massfunction(m,zz,rbs)

def dvdzdndmint(zz,Mmin,Mmax,omega_matter,h0,interp_da,rbs):
	Mint = scipy.integrate.quad(dndmint,np.log(Mmin),np.log(Mmax),args=(zz,rbs),epsrel=1e-6,epsabs = 0)[0]	
	return Mint*dVcdz(zz,omega_matter,h0,interp_da)

def dVcdz(z,omega_matter,h0,interp_da):
	c_light = 3.e5 #km/s
	da_z= interp_da(z)
	return c_light/(h0*100.0)*(da_z**2)*(1.0+z)**2/(np.sqrt(omega_matter*(1.0+z)**3 + (1.0 - omega_matter)))*(h0**3) # (Mpc/h)^3

def execute(block, config):
	# Configuration data, read from ini file above
	feedback,redshift,frac,minput,output_pdf  = config
	zmin = redshift - 0.01
	zmax = redshift + 0.01

	maxmass = block[cosmo, 'maxmass']
	omega_matter = block[cosmo, 'omega_m']
	h0 = block[cosmo, 'h0']

	z_da = block[dist,"z"][::-1]
	da_array = block[dist,"d_a"][::-1]
	interp_da = UnivariateSpline(z_da,da_array)
	zarray=block[mf,"z"]
	rarray=block[mf,"R_H"]
	rho_m=2.775e11*(omega_matter) # h^2 M_solar Mpc^-3.
	marray=(4.0*3.1415/3.0)*rho_m*rarray**3  # M_solar/h
	dndmarray=block[mf,"dndlnMh"].reshape([np.size(zarray),np.size(marray)]).T


	Mmin = 1.e12 #marray.min()
	Mmax = 1.e18 #marray.max()
	rbs = RectBivariateSpline(marray,zarray,dndmarray)

	ntot = scipy.integrate.quad(dvdzdndmint,zmin,zmax,args=(Mmin,Mmax,omega_matter,h0,interp_da,rbs),epsrel=1e-6,epsabs = 0)[0]
	NUM = ntot*frac

	
	LogPhi = np.zeros(minput.size)
	i=0
	if output_pdf:
		for mm in minput:
			FFm = frac/ntot*scipy.integrate.quad(dvdzdndmint,zmin,zmax,args=(Mmin,mm,omega_matter,h0,interp_da,rbs),epsrel=1e-6,epsabs = 0)[0]
			fm = frac/ntot*(scipy.integrate.quad(dvdm_zint,zmin,zmax,args=(mm,omega_matter,h0,interp_da,rbs),epsrel=1e-6,epsabs = 0)[0])
			LogPhi[i] = np.log(NUM*fm) + (NUM-1)*np.log(FFm)
			i = i + 1
		block[evs, 'logphi'] = LogPhi
		block[evs, 'm'] = minput


	FFm = frac/ntot*scipy.integrate.quad(dvdzdndmint,zmin,zmax,args=(Mmin,maxmass,omega_matter,h0,interp_da,rbs),epsrel=1e-6,epsabs = 0)[0]
	fm = frac/ntot*(scipy.integrate.quad(dvdm_zint,zmin,zmax,args=(maxmass,omega_matter,h0,interp_da,rbs),epsrel=1e-6,epsabs = 0)[0])
	LogLike = np.log(NUM*fm) + (NUM-1)*np.log(FFm)
	block[likes, 'EVS_LIKE'] = LogLike

		#signal that everything went fine
	return 0

def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness.  The joy of python.
    return 0
