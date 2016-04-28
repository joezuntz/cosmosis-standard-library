from cosmosis.datablock import option_section, names
import cl_likelihood
import numpy as np
import pdb

#TODO
## TO INCLUDE CMB KAPPA
#- Add kappa as additional observable with 1 tomographic bin
#- Read and interpolate noise from Tommaso to each ell
#- Read 1 extra survey: area, ell bins
#- Add override option to use given area rather than the smallest survey



def setup(options):

	import cl_likelihood as cll
	like = cll.ClLikelihood(options)

	cuts = options.get_bool(option_section, 'scale_cuts', default=False)

	shear_sample = options.get_string(option_section, 'shear_sample', default="")
	pos_sample = options.get_string(option_section, 'LSS_sample', default="")
	cmb_sample = options.get_string(option_section, 'cmb_sample', default="")

	auto = options.get_bool(option_section, 'auto_zbins')
	cross = options.get_bool(option_section, 'cross_zbins')

	override_area = options.get_bool(option_section, 'override_area', default=False)
	if override_area:
		area = options.get_double(option_section, 'area')

	save_dir = options.get_string(option_section, 'output', default="")

	config = [shear_sample, pos_sample, cuts], like, save_dir

	return config

def execute(block, config):

	options, like, output = config
	shear_cat = options[0]
	pos_cat = options[1]
#	cmb_cat = options[2]
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

	# Do the likelihood calculation
	like.do_likelihood(block)
	#like.normalise_likelihood(block)

	if output:
		np.savetxt(output, like.cov)
		print 'Saving covariance matrix to %s'%output

	# That's everything. The Cl likelihood class handles the details of
	# the covariance and likelihood calculations internally.
	return 0

def cleanup(config):
	return 0
