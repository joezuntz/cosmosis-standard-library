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

	#Specify gal_intrinsic_power=T to save gal_intrinsic power spectrum
	gal_intrinsic_power = options.get_bool(option_section,"gal_intrinsic_power",False)
	return method,gal_intrinsic_power

def execute(block, config):
	# load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
	lin = names.matter_power_lin
	nl = names.matter_power_nl
	ia = names.intrinsic_alignment_parameters
	gm = "galaxy_matter_power"
	cosmo = names.cosmological_parameters

	method,gal_intrinsic_power = config
	z_lin,k_lin,p_lin=block.get_grid(lin,"z","k_h","p_k")
	z_nl,k_nl,p_nl=block.get_grid(nl,"z","k_h","p_k")

	omega_m = block[cosmo, "omega_m"]
	A = block[ia, "A"]
	if abs(A)<1.e-9:
		A=1.e-9
	#if abs(A)<1e-6:
	#	block.put_grid(ia,"z",z_lin,"k_h",k_lin,"P_GI",np.zeros_like(p_nl))
	#	block.put_grid(ia,"z",z_lin,"k_h",k_lin,"P_II",np.zeros_like(p_nl))
	#run computation
	if method=='krhb':
		P_II, P_GI = kirk_rassat_host_bridle_power(z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)
	elif method=='bk':
		P_II, P_GI = bridle_king(z_nl, k_nl, p_nl, A, omega_m)
	elif method=='bk_corrected':
		P_II, P_GI = bridle_king_corrected(z_nl, k_nl, p_nl, A, omega_m)

	block.put_grid(ia, "z", z_nl, "k_h", k_nl, "p_matter_int", P_GI)
	block.put_grid(ia, "z", z_nl, "k_h", k_nl, "p_int_int", P_II)

	if gal_intrinsic_power:
		#This is a bit of hack...scale GI power spectrum (which is really matter-intrinsic
		#power spectrum by P_gal_matter/P_delta_delta
		z,k,p_gm=block.get_grid(gm, "z", "k_h", "p_k")
		P_gI = P_GI * p_gm/p_nl
		block.put_grid(ia, "z", z, "k_h", k, "p_pos_int", P_gI)
	return 0

def cleanup(config):
	pass
