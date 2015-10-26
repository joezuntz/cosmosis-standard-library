from cosmosis.datablock import names, option_section
from linear_alignments import kirk_rassat_host_bridle_power
from linear_alignments import bridle_king
from linear_alignments import bridle_king_corrected
import numpy as np

def setup(options):
	method = options[option_section, "method"].lower()
	if method not in ["krhb", "bk", "bk_corrected"]:
		raise ValueError('The method in the linear alignment module must'
			'be either "KRHB" (for Kirk, Rassat, Host, Bridle) or BK for '
			'Bridle & King or "BK_corrected" for the corrected version of that')
	return method

def execute(block, config):
	# load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
	lin = names.matter_power_lin
	nl = names.matter_power_nl
	ia = names.intrinsic_alignment_parameters
	cosmo = names.cosmological_parameters

	method = config

	z_lin,k_lin,p_lin=block.get_grid(lin,"z","k_h","p_k")
	z_nl,k_nl,p_nl=block.get_grid(nl,"z","k_h","p_k")

	omega_m = block[cosmo, "omega_m"]
	A = block[ia, "A"]

	#run computation and write to datablock
	if method=='krhb':
		P_II, P_GI, b_I, r_I = kirk_rassat_host_bridle_power(z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl,  "b_I", b_I)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl, "r_I", r_I)
	elif method=='bk':
		P_II, P_GI, b_I, r_I = bridle_king(z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia,  "z", z_nl, "k_h", k_nl, "b_I", b_I)
		block.put_grid(ia,  "z", z_nl, "k_h", k_nl, "r_I", r_I)
	elif method=='bk_corrected':
		P_II, P_GI, b_I, r_I = bridle_king_corrected(z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl, "b_I", b_I)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl,  "r_I", r_I)
		
	return 0

def cleanup(config):
	pass
