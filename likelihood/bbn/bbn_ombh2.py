from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

#The names of sections to read things from
cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


#Particle Data Group 2013
BBN_OMBH2_MEAN = 0.023
BBN_OMBH2_SIGMA = 0.002

def setup(options):
	mean = options.get_double(option_section, "mean", default=BBN_OMBH2_MEAN)
	sigma = options.get_double(option_section, "sigma", default=BBN_OMBH2_MEAN)
	norm = 0.5*log(2*pi*sigma**2)
	return (mean, sigma, norm)


def execute(block, config):
	# Configuration data, read from ini file above
	mean,sigma,norm = config


	# Get parameters from sampler
	h0 = block[cosmo, 'h0']
	omega_b = block[cosmo, 'omega_b']

	#compute the likelihood - just a simple Gaussian
	ombh2 = omega_b * h0**2
	like = -(ombh2-mean)**2/sigma**2/2.0 - norm
	block[likes, 'BBN_LIKE'] = like

	#signal that everything went fine
	return 0

def cleanup(config):
	#nothing to do here!  We just include this 
	# for completeness
	return 0
