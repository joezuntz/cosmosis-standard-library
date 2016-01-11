from cosmosis.datablock import names, option_section
import generate_observable_spectra as functions

import numpy as np
import matplotlib.pyplot as plt

def setup(options):
	shear= options[option_section, "shear"]
	intrinsic_alignments= options[option_section, "intrinsic_alignments"]
	clustering= options[option_section, "clustering"]
	magnification= options[option_section, "magnification"]
	noise = options[option_section, "noise"]
	bias = (options[option_section, "bias"], options[option_section, "m_per_bin"])
	bins = options[option_section, "angular_frequency_bins"]
	if bins:
		window = options[option_section, "window"]
	else:
		window = ''
	nlbin_shear = options[option_section, "nlbin_shear"]
	nlbin_ggl = options[option_section, "nlbin_ggl"]
	nlbin_pos = options[option_section, "nlbin_pos"]

	lmax_sh = options[option_section, "lmax_shear"]
	lmin_sh = options[option_section, "lmin_shear"]
	lmax_pos = options[option_section, "lmax_pos"]
	lmin_pos = options[option_section, "lmin_pos"]
	lmax_ggl = options[option_section, "lmax_ggl"]
	lmin_ggl = options[option_section, "lmin_ggl"]

	disp={True:'yes', False:'no'}
	print 'Shot noise: %s'%disp[noise]
	print 'Shape measurement bias (in each bin): %s (%s)'%(disp[bias[0]], disp[bias[1]])
	print 'Angular frequency binning: %s'%disp[bins]
	if bins:
		print 'Using window function: %s'%window

	shear_survey = options[option_section, "shear_sample"]
	pos_survey = options[option_section, "clustering_sample"]

	try: output_datavector = options[option_section, "output"]
	except: output_datavector = None

	opt= {'shear': shear, 'nlbin_shear': nlbin_shear, 'nlbin_ggl': nlbin_ggl, 'nlbin_pos': nlbin_pos, 'lmax_shear':lmax_sh, 'lmin_shear':lmin_sh, 'lmax_pos': lmax_pos, 'lmin_pos': lmin_pos, 'lmax_ggl': lmax_ggl, 'lmin_ggl': lmin_ggl, 'intrinsic_alignments': intrinsic_alignments, 'clustering': clustering, 'magnification': magnification, 'noise': noise, 'bias': bias, 'binning': bins, "window":window,'shear_cat': shear_survey, 'pos_cat': pos_survey, 'output_datavector': output_datavector}
	return opt

def execute(block, config):

	out_path = config['output_datavector']
	
	Cl = functions.Cl_class(block, config)
	Cl.load_and_generate_observable_cls(block, names)
	Cl.save_cls(block, out_path)

	return 0

def cleanup(config):
	pass
