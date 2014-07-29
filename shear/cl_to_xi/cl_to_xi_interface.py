#coding: utf-8
import cl_to_xi
import numpy as np
from cosmosis.datablock import option_section, names as section_names


def setup(options):
	config = {}
	try:
		n_theta = options[option_section, "n_theta"]
		theta_min = options[option_section, "theta_min"]
		theta_max = options[option_section, "theta_max"]
		theta_min = cl_to_xi.arcmin_to_radians(theta_min)
		theta_max = cl_to_xi.arcmin_to_radians(theta_max)
		theta = np.logspace(np.log10(theta_min), np.log10(theta_max), n_theta)
	except KeyError:
		try:
			theta = options[option_section, 'theta']
			if np.isscalar(theta): theta = [theta]
		except KeyError:
			print 'invalid theta input in options file for cl_to_xi'
			print 'need either theta = theta_1   theta_2  ... theta_n'
			print 'or theta_min, theta_max and n_theta.'
			print 'thetas in arcminutes'
			raise
	config['theta'] = theta
	return config


def execute(block, config):
	thetas = config['theta']
	n_theta = len(thetas)

	section = section_names.shear_cl
	output_section = section_names.shear_xi
	ell = block[section, "ell"]
	nbin = block[section, 'nbin']
	block[output_section, "theta"] = thetas
	block.put_metadata(output_section, "theta", "unit", "radians")

	for i in xrange(1,nbin+1):
		for j in xrange(1,i+1):
			name = 'bin_%d_%d'%(i,j)
			c_ell = block[section, name]
			xi_plus, xi_minus = cl_to_xi.cl_to_xi_plus_and_minus_extended_apodized(ell, c_ell, thetas)
			block[output_section, "xiplus_%d_%d"%(i,j)] = xi_plus
			block[output_section, "ximinus_%d_%d"%(i,j)] = xi_minus
	return 0

def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness.  The joy of python.
    return 0
