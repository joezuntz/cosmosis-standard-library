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

	config = [covmat_file, data_file, survey]

	return config

def execute(block, config):

	covmat_file = config[0]
	data_file = config[1]
	survey = config[2]

	pdb.set_trace()
	import cl_likelihood as cll
	like= cll.ClLikelihood(block, survey, covmat_file, data_file, shear=True, galaxy_clustering=True, sh_gal_cross=True, auto=True, cross=True)
	like.do_likelihood(block)
	pdb.set_trace()

	# That's everthing. The Cl likelihood class handles the details of
	# the covariance and likelihood calculations internally.
	return 0

def cleanup(config):
	return 0
