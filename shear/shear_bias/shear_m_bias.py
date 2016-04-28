"""

Errors in cosmic shear measurement can lead to a multiplicative factor
scaling the observed shear spectra.

This module scales the measured C_ell to account for that difference,
assuming model values of the multplicative factor m, either per bin or for all bins.

"""
from cosmosis.datablock import names, option_section

cal_section = names.shear_calibration_parameters

def setup(options):
	#This is an option - can set m_per_bin = T to get
	#a different m for each tomographic bin, or F to get
	#one global value
	m_per_bin=options.get_bool(option_section,"m_per_bin",True)
	return m_per_bin

def execute(block, config):

	m_per_bin=config
	if not m_per_bin:
		m0=block[cal_section, "m0"]

	cl_sec=names.shear_cl
	n_a=block[cl_sec,"nbin_a"]
	n_b=block[cl_sec,"nbin_b"]

	#Loop through bin pairs
	for i in xrange(1,n_a+1):
		for j in xrange(i,n_b+1):

			#Get existing C_ell
			cl_name="bin_%d_%d"%(j,i)
			cl_orig=block.get_double_array_1d(cl_sec,cl_name)

			#Compute scaling parameter on this pair
			if m_per_bin:
				mi = block[cal_section, "m%d"%i]
				mj = block[cal_section, "m%d"%j]
				m2 = (1+mi)*(1+mj)
			else:
				m2 = (1+m0)**2

			#Apply scaling and save back to block
			cl_new = m2*cl_orig
			block.replace_double_array_1d(cl_sec,cl_name,cl_new)

	return 0

	
