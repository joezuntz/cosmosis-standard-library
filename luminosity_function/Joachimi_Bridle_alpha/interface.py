from cosmosis.datablock import names, option_section
from luminosity_function import jb_calculate_alpha
import luminosity_function as luminosity
import numpy as np

def setup(options):

	sample = options.get_string(option_section, "sample", default='')
	binned_alpha = options.get_bool(option_section, "binned_alpha", default=True)

	if sample:
		mag_lim = options[sample, "magnitude_limit"]
	else:
		mag_lim = options[option_section, "magnitude_limit"]

	return (sample, binned_alpha, mag_lim )

def execute(block, config):

	ia = names.intrinsic_alignment_parameters
	cos = names.cosmological_parameters
	lum = names.galaxy_luminosity_function 
	
	survey, use_binned_alpha, mag_lim = config 

	# Obtain the fine grid points and limits for the redshift from the datablock
	# If the information isn't there already, set these parameters to sensible values
	try: Nz = int(block[survey, 'nz']) ;
	except: Nz = 500
	nzbin = int(block[survey, 'nzbin'])
	z = block[survey, 'z']		 

	coeff_a = luminosity.initialise_jb_coefficients(mag_lim)
	alpha, z = jb_calculate_alpha(coeff_a, z)

	# Write the fine grid alpha to the datablock
	try: block.put_double_array_1d(lum,'alpha',alpha)
	except IOError:
		print "Could not write alpha to datablock."
	block.put_double_array_1d(lum,'z',z)	

	# If required then interpolate alpha(z,rlim) to the mean values in each of the specified redshift bins
	if use_binned_alpha:
		alpha_bin, z_bar = luminosity.get_binned_alpha(block, survey, alpha, z)
	
		# Then write these to the datablock
                if len(alpha_bin)>1:
			block.put_double_array_1d(survey, 'alpha_binned', alpha_bin)
		else:
			block.put_double(survey, 'alpha_binned', alpha_bin[0])
		
	return 0

def cleanup(config):
	pass
