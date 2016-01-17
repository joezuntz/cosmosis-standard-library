from cosmosis.datablock import option_section, names
import scale_cuts
import numpy as np
import scipy.interpolate as inter
import pdb

def setup(options):
	allowed_cut_options = ['fixed', 'rassat08', 'file']

	scale_cut1 = options.get_string(option_section, 'method_shear').lower()
	scale_cut2 = options.get_string(option_section, 'method_LSS').lower()
	if scale_cut1 not in allowed_cut_options or scale_cut2 not in allowed_cut_options:
		raise ValueError("one of both of the scale cut options is recognised. Please choose from: 'fixed', 'Rassat08'")
		exit()

	if scale_cut1=='file' or scale_cut2=='file':
		filename = options.get_string(option_section, 'cuts_file')
	else:
		filename=None	

	cut_per_bin = None
	cut_per_bin = options.get_bool(option_section, 'cut_per_bin')

	shear = options.get_bool(option_section,'shear')
	ggl = options.get_bool(option_section,'ggl')
	pos = options.get_bool(option_section,'clustering')
	correlations = {'ee':shear, 'ne':ggl, 'nn':pos }

	shear_cat = options.get_string(option_section,'shear_sample')
	pos_cat = options.get_string(option_section,'clustering_sample')
	samples = {'shear': shear_cat, 'pos': pos_cat }

	return  correlations, samples, scale_cut1, scale_cut2, cut_per_bin, filename

def execute(block, config):
	corr,samples, cut_option_sh, cut_option_pos, cut_per_bin, filename = config

	# Set up a comoving distance interpolator if required
	chi_of_z_interpolator=None
	if cut_option_pos=='rassat08' or cut_option_sh=='rassat08':
		x = block['distances', 'd_m']
		z = block['distances', 'z']
		chi_of_z_interpolator = inter.interp1d(z,x)
		
	lmin_ee,lmax_ee,lmin_nn,lmax_nn, lmin_ne,lmax_ne = scale_cuts.get_angular_frequency_cuts(block, samples, cut_option_sh, cut_option_pos, cut_per_bin, corr, chi_of_z_interpolator, filename)
	scale_cuts.apply_angular_frequency_cuts(block, samples, corr, lmin_ee, lmax_ee, lmin_nn, lmax_nn, lmin_ne,lmax_ne)
	
	return 0

def cleanup(config):
	return 0
