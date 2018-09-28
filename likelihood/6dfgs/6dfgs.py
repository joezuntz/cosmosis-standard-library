from numpy import log, pi, interp, where, loadtxt,dot
import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
growthparams  = section_names.growth_parameters 
dist = section_names.distances


c_km_per_s =  299792.458


def setup(options):
	section = option_section
	mode = options.get_int(section, "mode", default=0)
	feedback = options.get_int(section, "feedback", default=0)

	if not mode:
		print("BAO likelihood")
		MEAN = 457. 
		SIGMA = 27.
		REDSHIFT = 0.106
	else:
		print("fsigma8 likelihood")
		MEAN = 0.423 
		SIGMA = 0.055
		REDSHIFT = 0.067
	
	mean = options.get_double(section, "mean", default=MEAN)
	sigma = options.get_double(section, "sigma", default=SIGMA)
	redshift = options.get_double(section, "redshift", default=REDSHIFT)
	norm = 0.5*log(2*pi*sigma**2)
	return (mode,mean,sigma,norm,redshift,feedback)


def execute(block, config):

	mode,mean,sigma,norm,redshift,feedback = config

	if not mode:
		dist_z = block[dist, 'z']
		d_m = block[dist, 'd_m']
		h = block[dist, 'h']
		if dist_z[1] < dist_z[0]:
			dist_z = dist_z[::-1]
			d_m = d_m[::-1]
			h = h[::-1]

		d_v = (dist_z*d_m*d_m/h)**(1./3.)
		d_v_predicted = interp(redshift, dist_z, d_v)
	
		if feedback:
			print("redshift = ", redshift)
			print("dv_predicted = ",  d_v_predicted)
			print("dv_data = ", mean)


		if (d_v_predicted<(mean+2.*sigma)) & (d_v_predicted>(mean-2.*sigma)):
			chi2 = (d_v_predicted-mean)**2/sigma**2
		else:
			chi2 = 4.

		like = -chi2/2.0 - norm

	else:
		z = block[growthparams, 'z']
		d_z = block[growthparams, 'd_z']
		f_z = block[growthparams, 'f_z']
		try:
			z0 = where(z==0)[0][0]
		except IndexError:
			raise ValueError("You need to calculate f(z) and d(z) down to z=0 to use the BOSS f*sigma8 likelihood")
		sig = block[cosmo, 'sigma_8']
		fsigma = (sig*(d_z/d_z[z0]))*f_z
		fsig = interp(redshift, z, fsigma)

		if feedback:
			print("Growth parameters: z = ",redshift, "fsigma_8  = ",fsig, " z0 = ", z0)

		like = -(fsig-mean)**2/sigma**2/2.0 - norm
	
	block[likes, '6DFGS_LIKE'] = like
	return 0
	

def cleanup(config):
    return 0
