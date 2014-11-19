from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

clusters = section_names.clusters
likes = section_names.likelihoods

MASS_MEAN = 4.02E14 # Santos et al. 2011
MASS_SIGMA = 0.88E14

def setup(options):
	section = option_section
	mean = options.get_double(section, "mean", default=MASS_MEAN)
	sigma = options.get_double(section, "sigma", default=MASS_SIGMA)
	norm = 0.5*log(2*pi*sigma**2)
	return (mean, sigma, norm)


def execute(block, config):
	# Configuration data, read from ini file above
	mean,sigma,norm = config

	# Get parameters from sampler
	maxmass = block[clusters, 'M_max']

	#compute the likelihood - just a simple Gaussian
	like = -(maxmass - mean)**2/sigma**2/2.0 - norm
	block[likes, 'MAXMASS_LIKE'] = like

	#signal that everything went fine
	return 0

def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness.  The joy of python.
    return 0
