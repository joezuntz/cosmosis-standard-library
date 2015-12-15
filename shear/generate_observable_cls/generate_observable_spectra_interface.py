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

	disp={True:'yes', False:'no'}
	print 'Shot noise: %s'%disp[noise]
	print 'Shape measurement bias (in each bin): %s (%s)'%(disp[bias[0]], disp[bias[1]])
	print 'Angular frequency binning: %s'%disp[bins]

	shear_survey = options[option_section, "shear_sample"]
	pos_survey = options[option_section, "clustering_sample"]

	try: output_datavector = options[option_section, "output"]
	except: output_datavector = None

	opt= {'shear': shear, 'intrinsic_alignments': intrinsic_alignments, 'clustering': clustering, 'magnification': magnification, 'noise': noise, 'bias': bias, 'binning': bins,'shear_cat': shear_survey, 'pos_cat': pos_survey, 'output_datavector': output_datavector}
	return opt

def execute(block, config):

	out_path = config['output_datavector']
	
	Cl = functions.Cl_class(block, config)
	Cl.load_and_generate_observable_cls(block, names)
	Cl.save_cls(block, out_path)

	return 0

def cleanup(config):
	pass
