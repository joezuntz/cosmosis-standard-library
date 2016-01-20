from cosmosis.datablock import names, option_section
from linear_alignments import kirk_rassat_host_bridle_power
from linear_alignments import bridle_king
from linear_alignments import bridle_king_corrected
import numpy as np

def setup(options):
	method = options[option_section, "method"].lower()
	spectra = options[option_section, "spectra"]
	colour ='none'
	try:
		colour = options.get_string(option_section, "colour").lower()
	except:
		colour='none'
	if method not in ["krhb", "bk", "bk_corrected"]:
		raise ValueError('The method in the linear alignment module must'
			'be either "KRHB" (for Kirk, Rassat, Host, Bridle) or BK for '
			'Bridle & King or "BK_corrected" for the corrected version of that')
	if colour not in ["none", "red", "blue", "r", "b"]:
		raise ValueError("Unrecognised choice of galaxy colour sample. Please specify"
				"'red' or 'blue' or omit the colour parameter in the config file.")
	return method, colour, spectra

def execute(block, config):
	# load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
	lin = names.matter_power_lin
	nl = names.matter_power_nl
	ia = names.intrinsic_alignment_parameters
	cosmo = names.cosmological_parameters

	method, colour, spectra = config

	z_lin,k_lin,p_lin=block.get_grid(lin,"z","k_h","p_k")
	z_nl,k_nl,p_nl=block.get_grid(nl,"z","k_h","p_k")

	omega_m = block[cosmo, "omega_m"]

	# Intended for backwards compatibility. If the user specifies a colour when the module
	# is called that is added to the name used in the datablock.
	# Otherwise we retain the old naming system and assume a single IA model for all galaxies
	if colour!='none':
		suffix = '_%s'%colour
	else:
		suffix = ''

	A = block[ia, "A"+suffix]

	#run computation and write to datablock
	if method=='krhb':
		P_II, P_GI, b_I, r_I, k_I = kirk_rassat_host_bridle_power(z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl,  "b_I"+suffix, b_I)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl, "r_I"+suffix, r_I)
		if spectra: 
			block.put_grid(ia, "z", z_lin, "k_h", k_I,  "P_GI"+suffix, P_GI)
			block.put_grid(ia, "z", z_lin, "k_h", k_I, "P_II"+suffix, P_II)
	elif method=='bk':
		P_II, P_GI, b_I, r_I = bridle_king(z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia,  "z", z_nl, "k_h", k_nl, "b_I"+suffix, b_I)
		block.put_grid(ia,  "z", z_nl, "k_h", k_nl, "r_I"+suffix, r_I)
		if spectra:
			block.put_grid(ia, "z", z_nl, "k_h", k_nl,  "P_GI"+suffix, P_GI)
                	block.put_grid(ia, "z", z_nl, "k_h", k_nl, "P_II"+suffix, P_II)
	elif method=='bk_corrected':
		P_II, P_GI, b_I, r_I = bridle_king_corrected(z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl, "b_I"+suffix, b_I)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl,  "r_I"+suffix, r_I)
		if spectra:
			block.put_grid(ia, "z", z_nl, "k_h", k_nl,  "P_GI"+suffix, P_GI)
                	block.put_grid(ia, "z", z_nl, "k_h", k_nl, "P_II"+suffix, P_II)

	return 0

def cleanup(config):
	pass
