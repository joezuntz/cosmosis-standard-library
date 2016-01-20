# Reads Schechter function fit parameters from a data file or gets
# them from the datablock
# Then calculates the luminosity function at the mean redshift of
# each tomographic bin
# Assumes that the bins are sufficiently narrow for the galaxy
# population in each to be characterised by a single set of Schechter 
# function parameters
# Simon Samuroff 12/15

from cosmosis.datablock import names, option_section
import numpy as np

valid_methods = ['faber_red', 'faber_blue', 'datablock']

def setup(options):
	survey = block[options, 'sample']
	mag = block[options, 'magnitude']
	lum = block[options, 'luminosity']
	method = block[options, 'method']
	per_bin = block[options, 'per_redshift_bin']

	config = survey, mag, lum, redshift, method, per_bin
	return config

def execute(block, config):

	sample, mag, lum, redshift, method, per_bin = config
	
	L = np.logspace(-25,1, 500)
	M = np.linspace(-25,1, 500)

	nzbin = block[sample, 'nzbin']

	if method=='faber_blue':
		data = np.loadtxt('data/faber_07_blue_fits.txt')
		z, alpha, phistar, Mstar=data.T
	elif method=='faber_red':
		data = np.loadtxt('data/faber_07_red_fits.txt')
		z, alpha, phistar, Mstar=data.T

	phi_of_M = {}
	if faber in method:
		z0 = block[sample, "z"]
		for bin in xrange(nzbin):
			if per_bin:
				# Get the mean redshift in this bin
				n_of_z = block[sample, "bin_%d"%(bin+1)]
				zf = np.linspace(z0[0], z0[-1], 1000)
				nz = np.interp(zf, z, n_of_z)
				z_mean = sum(nz*zf*(zf[1]-zf[0]))/sum(zf[1]-zf[0])

				# Interpolate to the mean redshift to find the Schechter function parameters
				a = np.interp(z_mean, z, alpha)
				phis = np.interp(z_mean, z, phistar)
				Ms = np.interp(z_mean, z, Mstar)
			else:
				a = alpha[0]
				phis = phistar[0]
				Ms = Mstar[0]

				phi_of_M['bin_'%(bin+1)] = 0.4 * np.log(10)* phis 
				phi_of_M['bin_'%(bin+1)] *= np.exp(-1.*10**(0.4*(Ms-M)) )
				phi_of_M['bin_'%(bin+1)] *= 10**(0.4*(Ms-M)*(a+1))
			

	elif method=='datablock':
		for bin in xrange(nzbin):
			if per_bin:
				a= block[sample, "alpha_%d"%(bin+1)]
				phis = block[sample, "phi_star_%d"%(bin+1)]
				if mag: Ms = block[sample, "M_star_%d"%(bin+1)]
				if lum: Ls = block[sample, "L_star_%d"%(bin+1)]
			else:
				a = block[sample, "alpha_1"]
				phis = block[sample, "phi_star_1"]
				if mag: Ms = block[sample, "M_star_1"]
				if lum: Ls = block[sample, "L_star_1"]
			if mag:
				phi_of_M['bin_'%(bin+1)] = 0.4 * np.log(10)* phis 
				phi_of_M['bin_'%(bin+1)] *= np.exp(-1.*10**(0.4*(Ms-M)) )
				phi_of_M['bin_'%(bin+1)] *= 10**(0.4*(Ms-M)*(a+1))
			if lum:
				phi_of_L['bin_'%(bin+1)] = phis * (L/Ls)**a  
				phi_of_L['bin_'%(bin+1)] *= np.exp(-1.*(L/Ls))

	#Save to datablock
	if lum:
		block.put_double_array_1d(sample, "L", L)
	if mag:
		block.put_double_array_1d(sample, "M", M) 
	for bin in nzbin:
		if lum:
			block[sample, "phi_of_L_%d"%(bin+1)] = phi_of_L["bin_%d"%(bin+1)]
		if mag:
			block[sample, "phi_of_M_%d"%(bin+1)] = phi_of_M["bin_%d"%(bin+1)]
	return 0

def cleanup(config):
	pass
