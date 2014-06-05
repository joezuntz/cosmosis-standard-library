import scipy.interpolate
import numpy as np
"""
This module generates the effective power P_12 
linear alignment KRHB model, which goes into C_ell
(under the Limber approximation) as :

C_ell = \int X_1(chi) X_2(chi) / chi^2 P_12(k=ell/chi,chi) d chi


"""
C1_BASELINE = 0.014
#in units of rho_crit0

def kirk_rassat_host_bridle_power(z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, A, Omega_m):
	""" 
	The Kirk, Rassat, Host, Bridle (2011) Linear Alignment model.
	Equations 8 and 9.

	C1 and rho_m0 must be in consistent units.
	The output P_II and P_GI will be specified at the same (k,z) as the linear power inputs

	The input f = A*C_1*rho_crit0

	"""

	#extrapolate our linear power out to high redshift
	z0 = np.where(z_lin==0)[0][0]
	nz = len(z_lin)

	# P_II is actually fixed across redshifts
	f = - Omega_m * A * C1_BASELINE

	# intrinsic-intrinsic term
	P_II = np.zeros_like(P_lin)
	for i in xrange(nz):
		P_II[:,i] = f**2 * P_lin[:,z0]

	#Get linear P(k) at the same sampling as non-linear
	P_nl_resample = np.zeros_like(P_lin)
	for i in xrange(nz):
		log_P_resample = np.interp(np.log(k_lin), np.log(k_nl), np.log(P_nl[:,i]))
		P_nl_resample[:,i] = np.exp(log_P_resample)

	growth = np.zeros_like(P_lin)
	ksmall = np.argmin(k_lin)
	for i in xrange(nz):
		growth[:,i] = (P_lin[ksmall,i] / P_lin[ksmall,z0])**0.5
	import pdb
	pdb.set_trace()
	P_GI = f * P_lin**0.5 * P_nl_resample**0.5 / growth

	return P_II, P_GI

