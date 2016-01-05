from cosmosis.datablock import option_section, names
from scipy.interpolate import interp1d
import numpy as np

MODES = ["multiplicative", "mean", "width"]

def setup(options):
	additive_bias = options[option_section, "mean"]
	broadening = options[option_section, "width"]
	survey = options[option_section, "survey"]
	if not additive_bias and not broadening and not multiplicative:
		raise ValueError("please set one or more of: %r to T"%MODES)
	return {"additive":additive_bias, "broadening":broadening, "survey":survey}

def execute(block, config):
	additive,broadening = config['additive'], config['broadening']
	pz = config['survey']
	biases = pz
	nbin = block[pz, "nzbin"]
	z = block[pz, "z"]
	for i in xrange(1,nbin+1):
		bin_name = "bin_%d" % i
		nz = block[pz, bin_name]
		bias = block[biases, "bias_%d"%i]
		dz = np.zeros_like(z)
		f = interp1d(z, nz, kind='cubic', fill_value = 0.0, bounds_error=False)
		#if mode=="multiplicative":
		#	nz_biased = f(z*(1-bias))
		if broadening:
			delta_z = block[biases, "delta_z_%d"%i]
			# Use the main peak of n(z) as a pivot point about which to distort n(z)
			zp = z[np.argwhere(nz==nz.max())[0][0]]
			dz += delta_z*(z-zp)
		if additive:
			dz -= bias
		nz_biased = f(z+dz)

		# Do a simple extrapolation to predict the points where n(z) has been set to zero
		# by the shift in the peak z or a change in width
		# Otherwise the n(z) is truncated at what was the end of the z array
		# This will start to look odd if the additive bias is large
		if nz_biased[-1]==0.0: 
			i = np.argwhere(nz_biased[0.5*len(nz_biased):]==0.0).T[0,0]
			alpha = (nz_biased[i-1]-nz_biased[i-2])/(z[i-1]-z[i-2])
			for j in xrange(len(z)-i):
				nzextrap = nz_biased[i+j-1] + alpha*(z[1]-z[0])
				if nzextrap>0:
					nz_biased[i+j] = nzextrap
				else:
					nz_biased[i+j] = 0.
		if nz_biased[0]==0.0: 
			i = np.argwhere(nz_biased[:0.5*len(nz_biased)]==0.0).T[0,-1]
			alpha = (nz_biased[i+1]-nz_biased[i+2])/(z[i+1]-z[i+2])
			for j in np.flipud(xrange(i)):
				nz_biased[j] = nz_biased[j+1] + alpha*(z[1]-z[0])
				nzextrap = nz_biased[j+1] + alpha*(z[1]-z[0])
				if nzextrap>0:
					nz_biased[j] = nzextrap
				else:
					nz_biased[j] = 0.

		#normalise
		nz_biased/=np.trapz(nz_biased,z)
		block[pz, bin_name] = nz_biased
	return 0

def cleanup(config):
	pass
