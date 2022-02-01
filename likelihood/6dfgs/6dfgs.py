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
	if options.has_value(section, "mode"):
		raise ValueError("The 6dfgs likelihood has changed and now uses parameters: "
						 "bao_likelihood, bao_choice, rsd_likelihood.  See the docs "
						 "for details")

	bao_like = options.get_bool(section, "bao_likelihood",default=True)
	bao_mode = options.get_string(section, "bao_choice",default='rs_dv')
	rsd_like = options.get_bool(section, "rsd_likelihood",default=False)
	
	if bao_like and rsd_like:
		raise ValueError("The RSD and BAO likelihoods can't be used together as we don't have their covariance matrix.")
	
	elif bao_like:
		if bao_mode == "dv":
			print("BAO likelihood on D_v")
			MEAN = 457. 
			SIGMA = 27.
			REDSHIFT = 0.106
		elif bao_mode == "rs_dv":
			print("BAO likelihood on r_s D_v")
			MEAN = 0.336 
			SIGMA = 0.015
			REDSHIFT = 0.106
		else:
			raise ValueError("Only the 6dfgs BAO likelihoods on D_v (keyword: 'dv') and r_s/D_v (keyword: 'rs_dv') are available")

	elif rsd_like:
		print("fsigma8 likelihood")
		MEAN = 0.423 
		SIGMA = 0.055
		REDSHIFT = 0.067
	
	else:
		raise ValueError("Specify the likelihood you want to use - BAO (bao_likelihood) or RSD (rsd_likelihood)")
	
	mean = options.get_double(section, "mean", default=MEAN)
	sigma = options.get_double(section, "sigma", default=SIGMA)
	redshift = options.get_double(section, "redshift", default=REDSHIFT)
	feedback = options.get_bool(section, "feedback", default=False)

	norm = 0.5*log(2*pi*sigma**2)
	return bao_like, bao_mode, mean, sigma, norm, redshift, feedback


def execute(block, config):

	bao_like, bao_mode, mean, sigma, norm, redshift, feedback = config

	if bao_like:
		dist_z = block[dist, 'z']
		d_m = block[dist, 'd_m']
		h = block[dist, 'h']
		if dist_z[1] < dist_z[0]:
			dist_z = dist_z[::-1]
			d_m = d_m[::-1]
			h = h[::-1]
		if bao_mode == 'rs_dv':
			rs_predicted = block[dist, "RS_ZDRAG"]
			z_drag_predicted = block[dist, "ZDRAG"]

		d_v = (dist_z*d_m*d_m/h)**(1./3.)
		d_v_predicted = interp(redshift, dist_z, d_v)
		
		if feedback:
			print("redshift = ", redshift)
			print("dv_predicted = ",  d_v_predicted)
			print("dv_data = ", mean)

		if bao_mode == 'dv':
			chi2 = (d_v_predicted - mean)**2 / sigma**2
		if bao_mode == 'rs_dv':
			rs_dv_predicted = rs_predicted / d_v_predicted
			chi2 = (rs_dv_predicted - mean)**2 / sigma**2


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
