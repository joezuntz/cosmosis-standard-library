from cosmosis.datablock import names, option_section
from suyu_likelihood import TimeDelayLikelihood, B1608, RXJ1131
import numpy as np

def setup(options):

	lens_name = options.get_string(option_section, "lens_name", "None")
	if lens_name.upper() == "B1608":
		like = B1608()
	elif lens_name.upper() == "RXJ1131":
		like = RXJ1131()
	elif lens_name=="None":
		try:
			z_d = options[option_section, "z_d"]
			z_s = options[option_section, "z_s"]
			lambda_d = options[option_section, "lambda_d"]
			mu_d = options[option_section, "mu_d"]
			sigma_d = options[option_section, "sigma_d"]
			like = TimeDelayLikelihood(z_d, z_s, lambda_d, mu_d, sigma_d)
		except:
			raise ValueError("If lens_name is not B1608 or RXJ1131, please set all of z_d, z_s, lambda_d, mu_d and sigma_d")

	return like


def execute(block, config):
	likelihood=config

	z_m = block[names.distances, "z"][::-1]
	d_m = block[names.distances, "d_m"][::-1]
	omega_k = block[names.cosmological_parameters, "omega_k"]
	H0 = block[names.cosmological_parameters, "hubble"]
	comovingDistance = lambda z: np.interp(z, z_m, d_m)

	like = likelihood.likelihood(comovingDistance, omega_k, H0)
	like_name = likelihood.name + "_LIKE"

	if np.isnan(like):
		like = -np.inf

	block[names.likelihoods, like_name] = like

	return 0



