"""

Errors in cosmic shear measurement can lead to a multiplicative factor
scaling the observed shear spectra.

This module scales the measured C_ell to account for that difference,
assuming model values of the multplicative factor m, either per bin or for all bins.

"""
from cosmosis.datablock import names, option_section

#cal_section = names.shear_calibration_parameters

def setup(options):
	#This is an option - can set m_per_bin = T to get
	#a different m for each tomographic bin, or F to get
	#one global value
	m_per_bin=options.get_bool(option_section,"m_per_bin",True)
	survey=options.get_string(option_section,"survey")
	return m_per_bin, survey

def execute(block, config):

	m_per_bin,survey=config
	cal_section = survey
	if not m_per_bin:
		m0=block[cal_section, "m0"]
		c0=block[cal_section, "c0"]

	cl_sec=names.shear_cl_gg
	n_z_bins=block[survey,"nzbin"]

	#Loop through bin pairs
	for i in xrange(1,n_z_bins+1):
		for j in xrange(i,n_z_bins+1):

			#Get existing C_ell
			cl_name="bin_%d_%d"%(j,i)
			cl_orig=block.get_double_array_1d(cl_sec,cl_name)

			#Compute scaling parameter on this pair
			if m_per_bin:
				mi = block[cal_section, "m%d"%i]
				mj = block[cal_section, "m%d"%j]
				m2 = (1+mi)*(1+mj)

				ci = block[cal_section, "c%d"%i]
				cj = block[cal_section, "c%d"%j]
				c2= ci*cj
			else:
				m2 = (1+m0)**2
				c2= c0*c0

			#Apply scaling and save back to block
			cl_new = m2*cl_orig + c2
			block.replace_double_array_1d(cl_sec,cl_name,cl_new)

	return 0

	
