from cosmosis.datablock import names
from cosmosis.datablock import option_section
import numpy as np
import scipy.integrate
from scipy.interpolate import RectBivariateSpline

#Option setup part.  Read options from in ifile.
# Definition of z-bins and R-bins
# Specify linear, nonlinear, or galaxy
def setup(options):
	if options.has_value(option_section, "z"):
		z = np.array(options[option_section, "z"])
	else:
		zmin = options[option_section, "zmin"]
		zmax = options[option_section, "zmax"]
		dz = options[option_section, "dz"]
		z = np.arange(zmin,zmax,dz)
	if options.has_value(option_section, "r"):
		R = np.array(options[option_section, "r"])
	else:
		rmin = options[option_section, "rmin"]
		rmax = options[option_section, "rmax"]
		dr = options[option_section, "dr"]		
		R = np.arange(rmin,rmax,dr)
	blockname = options[option_section, "matter_power"]
	print "Sigmar(R,z) will be evaluated at:"
	print "z = ", z
	print "R = ", R
	return (z, R, blockname)

def powerspec(k,z, rbs):
	return rbs.ev(k,z)
	#[0]
		
def sigint(lnk,r,z,rbs):
	k=np.exp(lnk)
	x=k*r
	w=3*(-x*np.cos(x)+np.sin(x))/x**3
	tmp = w**2*k**3*powerspec(k,z,rbs)/(2*3.14159**2)
	return tmp

def execute(block, config):
	z, R, blockname = config 

	zarray=block[blockname,"z"]
	karray=block[blockname,"k_h"]
	powerarray=block[blockname,"p_k"].reshape([np.size(zarray),np.size(karray)]).T

	rbs = RectBivariateSpline(karray,zarray,powerarray)

	kmin_overall = np.log(karray.min())
	kmax_overall = np.log(karray.max())

	sigma2r=np.zeros((np.size(R),np.size(z)))
	for i, rloop in enumerate(R):
		kmin=max(np.log(.01/rloop), kmin_overall)
		kmax=min(np.log(100./rloop), kmax_overall)
		for j, zloop in enumerate(z):
			sigma2r[i,j]=scipy.integrate.quad(sigint,kmin,kmax,args=(rloop,zloop,rbs),epsrel=1e-6)[0]
	section = "sigmar"	

	block.put_grid("sigma_r", "R", R, "z", z, "sigma2", sigma2r)
	return 0


def cleanup(config):
	#nothing to do here!  We just include this 
	# for completeness
	return 0
