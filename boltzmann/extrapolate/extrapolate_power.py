"""
This module takes linear and non-linear P(k) and extrapolates
them linearly in log-space out to a specified high k_max

"""
from cosmosis import option_section, names
import numpy as np
def extrapolate(log_k, log_P, n=200, log_new_kmax=None):
	""" 
	Extrapolate P(k) linearly in the logs of both P and k beyond
	where they are defined.  This is not specific to P(k) - it is
	just general linear interpolation.

	If log_new_kmax is set the the log_k must be the same for every redshift
	"""
	# Get the gradient at the end of the array
	log_k = log_k.squeeze()
	log_P = log_P.squeeze()
	nd = np.ndim(log_P)
	if nd==2:
		new_log_P = []
		d = log_P.shape[1]
		for i in xrange(d):
			log_P_i = log_P[:,i]
			new_log_k, new_log_P_i = extrapolate(log_k, log_P_i, n, log_new_kmax) 
			new_log_P.append(new_log_P_i)
		new_log_P = np.vstack(new_log_P).T
		return new_log_k, new_log_P

	elif nd>2:
		raise ValueError("Cannot")
	dlog_k = log_k[-1] - log_k[-2]
	dlog_P = log_P[-1] - log_P[-2]
	gradient = dlog_P / dlog_k

	if log_new_kmax is not None:
		n = int((log_new_kmax - log_k[-1])/dlog_k) + 1

	#Extend the k array with logarithmically spaced entries
	#and use the gradient we just found to extend P too.
	more_log_k = dlog_k * np.arange(n*1.0)
	more_log_P = more_log_k*gradient + log_P[-1]
	more_log_k += log_k[-1]
	#Combine the existing and extension arrays to get the
	#new extended array
	new_log_k = np.hstack((log_k, more_log_k))
	new_log_P = np.hstack((log_P, more_log_P))
	return new_log_k, new_log_P



def extrapolate_section(block, section, kmax):
	#load current values
	k = block[section, "k_h"]
	#early exit if we already have high enough k
	if k.max() >= kmax:
		return

	z = block[section, "z"]
	nk = len(k)
	nz = len(z)
	#load other current values
	P = block[section, "P_k"].reshape((nz, nk)).T
	#extrapolate
	logk, logp = extrapolate(np.log(k), np.log(P), log_new_kmax=np.log(kmax))
	#save results
	block[section, "k_h"] = np.exp(logk)
	block[section, "P_k"] = np.exp(logp).T.flatten()
	block[section, "nk"] = len(logk)


def setup(options):
	kmax = options.get_double(option_section, "kmax")
	return {"kmax":kmax}


def execute(block, config):
	kmax = config['kmax']
	#extrapolate non-linear power
	for section in [names.matter_power_nl, names.matter_power_lin]:
		if block.has_section(section):
			extrapolate_section(block, section, kmax)
	return 0
