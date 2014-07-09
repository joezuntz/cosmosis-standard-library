from cosmosis import option_section, names
import clerkin
import numpy as np

models = ['gtd', 'q', 'q-gtd'] #...]
modes = ["power","bias","both"]
def setup(options):
	model = options[option_section, "model"].lower()
	mode = options[option_section, "mode"].lower()
	if model not in models:
		raise ValueError("The Clerkin module requires a model chosen from: %r" % models)
	if mode not in modes:
		raise ValueError("The Clerkin module requires a mode, chosen from: %r" % modes)
	return model, mode

def execute_power(block, model):

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
def execute_bias(block, model):
	if model=="gtd":
		#Get growth
		z1 = block[names.growth_parameters,'z']
		growth_z = block[names.growth_parameters,'d_z']
		#Get bias parameters
		alpha = block['bias_parameters', 'alpha']
		b0 = block['bias_parameters', 'b0']
		c = block['bias_parameters', 'c']
		b = clerkin.gtd_bias(z1, growth_z, alpha, b0, c)
		block[names.bias_field, "z"] = z1
		block[names.bias_field, "b"] = b

	elif model=="q":
		k = np.logspace(-6, 2, 500)	
		Q = block['bias_parameters', 'Q']
		A = block['bias_parameters', 'A']
		b = clerkin.q_bias(k, Q, A)
		block[names.bias_field, "k"] = k
		block[names.bias_field, "b"] = b

	elif model=='q-gtd':
		z1 = block[names.growth_parameters,'z']
		growth_z = block[names.growth_parameters,'d_z']
		#Get bias parameters
		alpha = block['bias_parameters', 'alpha']
		b0 = block['bias_parameters', 'b0']
		c = block['bias_parameters', 'c']
		b = clerkin.gtd_bias(z1, growth_z, alpha, b0, c)
		k = np.logspace(-6, 2, 500)	
		Q = block['bias_parameters', 'Q']
		A = block['bias_parameters', 'A']
		b_z = clerkin.gtd_bias(z1, growth_z, alpha, b0, c)
		b_k = clerkin.q_bias(k, Q, A)
		nk = len(k)
		nz = len(z1)
		b = np.zeros((nk,nz))
		for i in xrange(nk):
			for j in xrange(nz):
				b[i, j] = b_k[i]*b_z[j]
		block.put_grid(names.bias_field, "k_h", k, "z", z1, "b", b)
	return 0

def execute_both(block, model):
	status = execute_bias(block, model)
	if status: return status
	return execute_power(block, model)

def execute(block, config):
	model,mode = config
	if mode=="bias":
		return execute_bias(block, model)
	elif mode=="power":
		return execute_power(block, model)
	elif mode=="both":
		return execute_both(block, model)
	else:
		raise ValueError("Unknown mode in Clerkin")


