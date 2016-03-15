from cosmosis.datablock import option_section, names
import cl_likelihood
import numpy as np
import pdb

def setup(options):

	import cl_likelihood as cll
	like = cll.ClLikelihood(options)

	cuts = options.get_bool(option_section, 'scale_cuts')


	shear_sample = options.get_string(option_section, 'shear_sample')
	pos_sample = options.get_string(option_section, 'LSS_sample')
	auto = options.get_bool(option_section, 'auto_zbins')
	cross = options.get_bool(option_section, 'cross_zbins')

	try: save_dir = options.get_string(option_section, 'output')
	except: save_dir = None

	config = [shear_sample, pos_sample, cuts], like, save_dir

	return config

def execute(block, config):

	options, like, output = config
	shear_cat = options[0]
	pos_cat = options[1]
	cuts = options[2]

	# First apply scale cuts
	# This is done here rather than on setup as they may change with
	# each iteration of the pipeline (if lmin, lmax are fitted as
	# free parameters)
	like.apply_scale_cuts(block, cuts)
	# Setup the theory vector
	like.initialise_theory(block)

	if like.constant_covariance:
		like.build_inverse_covariance(block)
	#import pdb ; pdb.set_trace()

	# Do the likelihood calculation
	like.do_likelihood(block)
	#like.normalise_likelihood(block)

	if output!=None:
		np.savetxt(output, like.cov)
		print 'Saving covariance matrix to %s'%output

	# That's everything. The Cl likelihood class handles the details of
	# the covariance and likelihood calculations internally.
	return 0

def cleanup(config):
	return 0
