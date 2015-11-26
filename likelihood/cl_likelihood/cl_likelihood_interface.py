from cosmosis.datablock import option_section, names
import cl_likelihood
import numpy as np
import pdb

def setup(options):
	sec = option_section
	try: covmat_file = options.get_string(sec, 'covariance')
	except: covmat_file=None
	data_file = options.get_string(sec, 'data')
	survey = options.get_string(sec, 'survey')
	shear = options.get_bool(sec, 'shear')
	galaxy_clustering = options.get_bool(sec, 'galaxy_clustering')
	ggl = options.get_bool(sec, 'ggl')
	auto = options.get_bool(sec, 'auto_zbins')
	cross = options.get_bool(sec, 'cross_zbins')

	import cl_likelihood as cll
	like= cll.ClLikelihood(survey, covmat_file, data_file, shear=shear, galaxy_clustering=galaxy_clustering, sh_gal_cross=ggl)

	config = [covmat_file, data_file, survey, auto, cross], like

	return config

def execute(block, config):

	options, like = config
	covmat_file = options[0]
	data_file = options[1]
	survey = options[2]
	auto = options[3]
	cross = options[4]
	like.initialise_theory(block, covmat_file, auto=auto, cross=cross)
	like.do_likelihood(block)

	# That's everthing. The Cl likelihood class handles the details of
	# the covariance and likelihood calculations internally.
	return 0

def cleanup(config):
	return 0
