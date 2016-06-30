from cosmosis.datablock import names, option_section
from time_delay_likelihood import TimeDelayLikelihood, B1608, RXJ1131
import numpy as np

def setup(options):

	lens_name = options.get_string(option_section, "lens_name", "B1608")
	if lens_name.upper() == "B1608":
		like = B1608()
	elif lens_name.upper() == "RXJ1131":
		like = RXJ1131()
	elif lens_name=="mock":
		try:
                        fname = options[option_section, "filename"]
                        lambdad = options[option_section, "lambdaD"]
                        #mock_data = np.loadtxt("fname",skiprows = 1)
			#z_d = mock_data[:,0]; z_s = mock_data[:,1] ; d_obs = mock_data[:,2]; sigma_d = mock_data[:,2]
			#lambda_d = np.zeros(len(d_obs))*4000.
			#mu_d = np.log(d_obs - lambda_d)
                        like = TimeDelayLikelihood.load_catalog(fname,lambdad) # rtns array of instances of class
		except:
			raise ValueError("Error in reading mocks: z_d, z_s, lambda_d, mu_d and sigma_d")

	return like


def execute(block, config):
	data_class=config

	z_m = block[names.distances, "z"][::-1]
	d_m = block[names.distances, "d_m"][::-1]
	omega_k = block[names.cosmological_parameters, "omega_k"]
	H0 = block[names.cosmological_parameters, "hubble"]
	comovingDistance = lambda z: np.interp(z, z_m, d_m)

        if len(data_class)>1:
                like = 0 
                for d in data_class:
#                    print "LIKE", like
                    like = like + d.likelihood(comovingDistance, omega_k, H0)
                    like_name = d.name + "_LIKE"

        else:
	    like = data_class.likelihood(comovingDistance, omega_k, H0)
	    like_name = data_class.name + "_LIKE"

	if np.isnan(like):
		like = -np.inf

	block[names.likelihoods, like_name] = like

	return 0



