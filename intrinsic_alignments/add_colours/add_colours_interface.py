from cosmosis.datablock import option_section, names
import numpy as np

def setup(options):
	red_model = options[option_section, "red_model"]
	blue_model = options[option_section, "blue_model"]
	catalogue = options[option_section, "catalogue"]
	return red_model, blue_model, catalogue

def execute(block, config):
	modelr, modelb, catalogue = config

	section_II = names.ia_spectrum_ii
	section_GI = names.ia_spectrum_gi
	section_ia = names.intrinsic_alignment_parameters

	# Load nuisance parameters for the two models
	f_red = block[catalogue, 'red_fraction']
	f_blue = 1. - f_red
	A_red = block[section_ia, 'A' + '_%s'%model_r]
	A_blue = block[section_ia, 'A' + '_%s'%model_b]
	alpha_red = block[section_ia, 'alpha' + '_%s'%model_r]
	alpha_blue = block[section_ia, 'alpha' + '_%s'%model_b]

	# Define scaling grid
	_,z_grid=np.meshgrid(k,z)

	# Get relevant power spectra
	z,k,P_II_red = block.get_grid(section_ia,  "z", "k_h", "P_II"+ '_%s'%model_r)
	z,k,P_II_blue = block.get_grid(section_ia,  "z", "k_h", "P_II"  + '_%s'%model_b)
	z,k,P_GI_red  = block.get_grid(section_ia,  "z", "k_h", "P_GI" + '_%s'%model_r)
	z,k,P_GI_blue  = block.get_grid(section_ia,  "z", "k_h", "P_GI" + '_%s'%model_b)
	P_II_red_blue = np.sqrt(P_II_red*P_II_blue)

	# Combine red, blue and cross terms
	P_II = f_red*f_red * A_red*A_red * P_II_red * (1+z)**(2.*alpha_red)
	P_II += f_blue*f_blue * A_blue*A_blue * P_II_blue * (1+z)**(2.*alpha_blue)
	P_II += f_red*f_blue * A_red*A_blue * P_II_red_blue * (1+z)**(alpha_blue)*(1+z)**(alpha_red)

	P_GI = f_red * A_red * P_GI_red * (1+z)**(alpha_red)
	P_GI = f_blue * A_blue * P_GI_blue * (1+z)**(alpha_blue)
