from cosmosis import names, option_section
# sigma8 * omega_m ** 0.3 = 

cosmo = names.cosmological_parameters
likes = names.likelihoods
sz = "sz"

def setup(options):
	measured_value = options.get_double(option_section, "measured_value", 0.764)
	error = 0.025	
	fiducial_omega = 0.27
	return (measured_value, fiducial_omega, error)

def execute(block, config):
	(mu, fid, sigma) = config
	sigma8 = block[cosmo, "sigma_8"]
	omega_m = block[cosmo, "omega_m"]

	# Cosmology ...
	x = sigma8 * (omega_m/fid)**0.3

	# save the likelihood
	like = -0.5 * (x-mu)**2 / sigma**2
	block[likes, "SZ_LIKE"] = like
	block[sz, "x"] = x
	#could put:
	#extra_outputs = sz/x

	return 0


def cleanup(config):
	pass

