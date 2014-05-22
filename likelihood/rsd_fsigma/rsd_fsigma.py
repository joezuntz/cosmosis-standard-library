from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
growthparams  = section_names.growth_parameters 

fsig_MEAN = 0.428 # Chuang et al 2013 BOSS DR9
fsig_SIGMA = 0.066
REDSHIFT = 0.57

def setup(options):
	section = option_section
	mean = options.get_double(section, "mean", default=fsig_MEAN)
	sigma = options.get_double(section, "sigma", default=fsig_SIGMA)
	redshift = options.get_double(section, "redshift", default=REDSHIFT)
	feedback = options.get_int(section, "feedback", default=0)
	norm = 0.5*log(2*pi*sigma**2)
	return (mean, sigma, norm,redshift,feedback)


def execute(block, config):
	# Configuration data, read from ini file above
	mean,sigma,norm,redshift,feedback = config

	# Get parameters from sampler
	Dz = block[growthparams, 'D_z']
	Dz0 = block[growthparams, 'D_z0']
	fz = block[growthparams, 'f_z']
	z = block[growthparams, 'growth_z']
	sig = block[cosmo, 'sigma_8']
	if redshift != z:
		print "Error : redshift for growth module and  RSD fsig are different in ini file."
		return 1
	fsig = (sig*(Dz/Dz0)**2)*fz
	if feedback:
		print "Growth parameters: z = ",redshift, "fsigma_8  = ",fsig, "D = ",Dz, "f = ",fz
	#compute the likelihood - just a simple Gaussian
	like = -(fsig-mean)**2/sigma**2/2.0 - norm
	block[likes, 'RSD_FSIGMA_LIKE'] = like

	#signal that everything went fine
	return 0

def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness.  The joy of python.
    return 0
