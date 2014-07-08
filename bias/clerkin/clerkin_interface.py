from cosmosis import option_section, names
import clerkin

models = ['gtd', 'q', 'q-gtd'] #...]

def setup(options):
	model = options[option_section, "model"].lower()
	if model not in models:
		raise ValueError("The Clerkin module requires a model chosen from: %r" % models)
	return model

def execute(block, config):
	model = config


	#Get matter power
	k, z, P = block.get_grid(names.matter_power_nl, "k_h", "z", "p_k")



	#Get matter power with bias
	if model=='gtd':
		#Get growth
		z1 = block[names.growth_parameters,'z']
		growth_z = block[names.growth_parameters,'d_z']
		#Get bias parameters
		alpha = block['bias_parameters', 'alpha']
		b0 = block['bias_parameters', 'b0']
		c = block['bias_parameters', 'c']
		#compute new power
		P_out = clerkin.gtd_model(z1, growth_z, k, z, P, alpha, b0, c)
	elif model=='q':
		#Get bias parameters
		Q = block['bias_parameters', 'Q']
		A = block['bias_parameters', 'A']
		#compute new power
		P_out = clerkin.q_model(k, z, P, Q, A)
	elif model=='q-gtd':
		#Get growth
		z1 = block[names.growth_parameters,'z']
		growth_z = block[names.growth_parameters,'d_z']

		#Get bias parameters
		alpha = block['bias_parameters', 'alpha']
		b0 = block['bias_parameters', 'b0']
		c = block['bias_parameters', 'c']
		Q = block['bias_parameters', 'Q']
		A = block['bias_parameters', 'A']
		#compute power
		P_out = clerkin.gtd_q_model(z1, growth_z, k, z, P, alpha, b0, c, Q, A)
	else:
		raise ValueError("Unknown error in Clerkin")
	
	if (P_out<0).any():
		print "Negative P: ", P_out.min()
		return 1

	block.put_grid("matter_power_gal", "z", z, "k_h", k, "p_k", P_out.T)

	return 0