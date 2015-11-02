from cosmosis.datablock import option_section, names
from scipy.interpolate import interp1d
import numpy as np

MODES = ["multiplicative", "additive"]

def setup(options):
	mode = options[option_section, "mode"]
	if mode not in MODES:
		raise ValueError("mode for photoz must be one of: %r"%MODES)
	return {"mode":mode}

def execute(block, config):
	mode = config['mode']
	pz = names.wl_number_density
	biases = "wl_photoz_errors"
	nbin = block[pz, "nbin"]
	z = block[pz, "z"]
	for i in xrange(1,nbin+1):
		bin_name = "bin_%d" % i
		nz = block[pz, bin_name]
		bias = block[biases, "bias_%d"%i]
		f = interp1d(z, nz, kind='cubic', fill_value = 0.0, bounds_error=False)
		if mode=="multiplicative":
			nz_biased = f(z*(1-bias))
		elif mode=="additive":
			nz_biased = f(z-bias)
		else:
			raise ValueError("Unknown photo-z mode")
		#normalize
		nz_biased/=np.trapz(nz_biased,z)
		block[pz, bin_name] = nz_biased
	return 0

def cleanup(config):
	pass
