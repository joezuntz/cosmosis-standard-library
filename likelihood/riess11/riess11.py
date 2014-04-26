from numpy import log, pi
import cosmosis.datablock.cosmosis_py


cosmo = cosmosis_py.section_names.cosmological_parameters
likes = cosmosis_py.section_names.likelihoods

HST_H0_MEAN = 0.738
HST_H0_SIGMA = 0.024

def setup(options):
	section = cosmosis_py.block.option_section
	mean = options.get_double(section, "mean", default=HST_H0_MEAN)
	sigma = options.get_double(section, "sigma", default=HST_H0_SIGMA)
	norm = 0.5*log(2*pi*sigma**2)
	return (mean, sigma, norm)


def execute(block, config):
	# Configuration data, read from ini file above
	mean,sigma,norm = config

	# Get parameters from sampler
	h0 = block[cosmo, 'h0']

	#compute the likelihood - just a simple Gaussian
	like = -(h0-mean)**2/sigma**2/2.0 - norm
	block[likes, 'RIESS_LIKE'] = like

	#signal that everything went fine
	return 0
