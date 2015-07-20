#IA model from maccrann et al. 2015 - NLA with power law redshift scaling
#P_II(k,z) = ((1+z)^alpha) * F(z))^2 * P_delta(k,z), P_GI(k,z) = (1+z)^alpha * F(z) * P_delta(k,z)
#The P_II and P_GI before the redshift scaling are already in the block
#so just read them in and rescale

from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
	return 0

def execute(block, config):

	ia_section = names.intrinsic_alignment_parameters

	#read in power spectra
	z,k,p_ii=block.get_grid(ia_section,"z","k_h","P_II")
	z,k,p_gi=block.get_grid(ia_section,"z","k_h","P_GI")

	#read alpha from ia_section values section - JAZ this is a change from the internal one
	alpha = block[ia_section,'alpha']
	_,z_grid=np.meshgrid(k,z)

	#Construct and Apply redshift scaling
	z_scaling=(1+z_grid)**alpha
	p_ii*=z_scaling**2
	p_gi*=z_scaling

	#Save grid
	block.replace_grid(ia_section, "z", z, "k_h", k, "P_GI", p_gi)
	block.replace_grid(ia_section, "z", z, "k_h", k, "P_II", p_ii)
	return 0

def cleanup(config):
	pass