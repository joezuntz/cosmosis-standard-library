import numpy as np
import cosmosis_py.section_names
from cosmosis_py.block import option_section

def setup(options):
	#only one parameter - filepath
	filename = options[option_section, "filepath"]
	data = np.loadtxt(filename).T
	nz = len(data[0])
	nbin = len(data)-1
	z = data[0]
	n_of_z = data[1:]

	#Normalize n(z)
	for col in n_of_z:
		norm = np.trapz(col, z)
		col/=norm

	print "Found %d samples and %d bins in redshift in file %s" % (nbin, nz, filename)
	return (nz, nbin, z, n_of_z)

def execute(block, config):
	(nz, nbin, z, n_of_z) = config

	section = cosmosis_py.section_names.wl_number_density
	block[section, 'nz'] = nz
	block[section, 'nbin'] = nbin
	block[section, 'z'] = z
	for (bin, bin_n_of_z) in enumerate(n_of_z):
		name = "bin_%d"%(bin+1)
		block[section, name] =  bin_n_of_z

	return 0


def cleanup(config):
	#nothing to do here!  We just include this 
	# for completeness
	return 0
