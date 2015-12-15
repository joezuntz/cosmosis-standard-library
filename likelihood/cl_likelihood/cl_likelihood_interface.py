from cosmosis.datablock import option_section, names
import cl_likelihood
import numpy as np
import pdb

def setup(options):

	import cl_likelihood as cll
	like = cll.ClLikelihood(options)

	shear_sample = options.get_string(option_section, 'shear_sample')
	pos_sample = options.get_string(option_section, 'clustering_sample')
	auto = options.get_bool(option_section, 'auto_zbins')
	cross = options.get_bool(option_section, 'cross_zbins')

	config = [shear_sample, pos_sample, auto, cross], like

	return config

def execute(block, config):

	options, like = config
	shear_cat = options[0]
	pos_cat = options[1]
	auto = options[2]
	cross = options[3]

	# First apply scale cuts
	# This is done here rather than on setup as they may change with
	# each iteration of the pipeline (if lmin, lmax are fitted as
	# free parameters)
	like.apply_scale_cuts(block)
	# Then calculate the covariance matrix
	if like.constant_covariance:
		like.extract_covariance(block)
	# Setup the theory vector
	like.initialise_theory(block, auto=auto, cross=cross)
	# Do the likelihood calculation
	like.do_likelihood(block)

	# That's everything. The Cl likelihood class handles the details of
	# the covariance and likelihood calculations internally.
	return 0

def cleanup(config):
	return 0
