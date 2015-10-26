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

	survey = options[option_section, "survey"]

	try: output_datavector = options[option_section, "output"]
	except: output_datavector = None

	opt= {'shear': shear, 'intrinsic_alignments': intrinsic_alignments, 'clustering': clustering, 'magnification': magnification, 'noise': noise, 'survey': survey, 'output_datavector': output_datavector}
	return opt

def execute(block, config):
	
	survey = config['survey']

	# Survey parameters
	Nl = int(block[survey, 'nlbin'])
	lmax = block.get_double(survey, 'lmax')
	lmin = block.get_double(survey, 'lmin')

	# Defines whether to add noise 
	noise = config['noise']

	out_path = config['output_datavector']
	
	Cl = functions.Cl_class(block, config)
	Cl.load_and_generate_observable_cls(block, names)
	Cl.save_cls(block, out_path)

	return 0

def cleanup(config):
	pass
