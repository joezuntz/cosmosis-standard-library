from cosmosis.datablock import names, option_section
from linear_alignments import kirk_rassat_host_bridle_power
from linear_alignments import bridle_king

def setup(options):
	method = options[option_section, "method"].lower()
	if method not in ["krhb", "bk"]:
		raise ValueError('The method in the linear alignment module must'
			'be either "KRHB" (for Kirk, Rassat, Host, Bridle) or BK for '
			'Bridle & King')
	return method

def execute(block, config):
	# load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
	lin = names.matter_power_lin
	nl = names.matter_power_nl
	ia = names.intrinsic_alignment_parameters
	cosmo = names.cosmological_parameters

	method = config

	z_lin = block[lin, "z"]
	k_lin = block[lin, "k_h"]
	p_lin = block[lin, "p_k"]
	z_nl = block[nl, "z"]
	k_nl = block[nl, "k_h"]
	p_nl = block[nl, "p_k"]

	p_lin = p_lin.reshape((z_lin.size, k_lin.size)).T
	p_nl = p_nl.reshape((z_nl.size, k_nl.size)).T


	omega_m = block[cosmo, "omega_m"]
	A = block[ia, "A"]

	#run computation
	if method=='krhb':
		P_II, P_GI = kirk_rassat_host_bridle_power(z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia, "z", z_lin, "k_h", k_lin, "P_GI", P_GI.T)
		block.put_grid(ia, "z", z_lin, "k_h", k_lin, "P_II", P_II.T)
	elif method=='bk':
		P_II, P_GI = bridle_king(z_nl, k_nl, p_nl, A, omega_m)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl, "P_GI", P_GI.T)
		block.put_grid(ia, "z", z_nl, "k_h", k_nl, "P_II", P_II.T)



	return 0

def cleanup(config):
	pass