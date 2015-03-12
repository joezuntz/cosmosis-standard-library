import numpy as np
from cosmosis.datablock import option_section, names as section_names

def setup(options):
	#only one parameter - filepath
	filename = options[option_section, "filepath"]
	des_fmt = options.get_bool(option_section,"des_fmt",default=False)
	data_full = np.loadtxt(filename).T
	if des_fmt:
		z=0.5*(data_full[0]+data_full[1])
		nz=len(z)
		nbin=len(data_full)-3
		n_of_z=data_full[2:-1]
	else:
		nz = len(data_full[0])
		nbin = len(data_full)-1
		z = data_full[0]
		n_of_z = data_full[1:]
	#check first z is zero, if not add some
	if z[0]>0.00000001:
		z_new=np.zeros(len(z)+1)
		z_new[1:]=z
		n_of_z_new=np.zeros((nbin,len(z)+1))
		n_of_z_new[:,1:]=n_of_z
		z,n_of_z=z_new,n_of_z_new

	#Normalize n(z)
	for col in n_of_z:
		norm = np.trapz(col, z)
		col/=norm

	print "Found %d samples and %d bins in redshift in file %s" % (nbin, nz, filename)
	return (nz, nbin, z, n_of_z)

def execute(block, config):
	(nz, nbin, z, n_of_z) = config

	section = section_names.wl_number_density
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
