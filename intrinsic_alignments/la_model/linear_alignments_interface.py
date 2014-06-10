from cosmosis import names, option_section
from linear_alignments import kirk_rassat_host_bridle_power

def setup(options):
	return 1

def execute(block, config):
	# load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
	lin = names.matter_power_lin
	nl = names.matter_power_nl
	ia = names.intrinsic_alignment_parameters
	cosmo = names.cosmological_parameters

	II = "intrinsic_alignment_II"
	GI = "intrinsic_alignment_GI"

	z_lin = block[lin, "z"]
	k_lin = block[lin, "k_h"]
	p_lin = block[lin, "p_k"]
	z_nl = block[lin, "z"]
	k_nl = block[lin, "k_h"]
	p_nl = block[lin, "p_k"]

	p_lin = p_lin.reshape((z_lin.size, k_lin.size)).T
	p_nl = p_nl.reshape((z_nl.size, k_nl.size)).T

	# import pylab
	# pylab.loglog(k_lin, p_lin[:,0])
	# pylab.loglog(k_nl, p_nl[:,0])

	omega_m = block[cosmo, "omega_m"]
	A = block[ia, "A"]

	#run f computation
	P_II, P_GI = kirk_rassat_host_bridle_power(z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)

#	pylab.loglog(k_lin, P_II[:,0])
	pylab.show()

	block.put_grid(II, "k_h", k_lin, "z", z_lin, "P_k", P_II)
	block.put_grid(GI, "k_h", k_lin, "z", z_lin, "P_k", P_GI)

	return 0

def cleanup(config):
	pass