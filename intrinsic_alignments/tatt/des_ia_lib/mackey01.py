#Mackey, White, Kamionkowski 0106364 model
#Eqn 27 - we just need to get P_II - this will get limbered later by the spectra module
import scipy.interpolate as interp
import numpy as np
import scipy.integrate
import del4
#import pylab

pi=np.pi
"""
def intrinsic_spectra(z_lin, k_lins, P_lin, z_growth, Dz):

	p_lin_z0=P_lin[0]
	#We need the growth function to rescale the power spectrum at the end...get growth at same z-values as P_lin
	assert np.isclose(z_growth[0],0.)
	assert ((z_growth[-1]>z_lin[-1]) or np.isclose(z_growth[-1],z_lin[-1]))
	if not np.allclose(z_lin, z_growth):
		Dz=interp.interp1d(z_growth,Dz)(z_lin)

	alphas=np.sort(np.concatenate((np.logspace(-4,4,50),np.linspace(0.98,1.02,50))))
	mus=np.concatenate((np.linspace(-1.,0.99,50),np.linspace(0.99,1.,50)))

	P_ee = np.zeros_like(p_lin_z0)

	del2_interp = del4.del2_k_interp_func(k_lins, p_lin_z0)
	R=1
	sig2R=del2_interp.sigma2_r(R)
	#print 'sig2R',sig2R
	pref = 5 * 0.4**2  / 32 / (sig2R**2)

	smooth=1.
	alpha_integrand=np.zeros((len(alphas),len(k_lins)))
	for ialpha,alpha in enumerate(alphas):
		mu_integrand=np.zeros((len(mus),len(k_lins)))
		for imu,mu in enumerate(mus):
			polyee = g_ee(alpha,mu)
			k_diff_ang=np.sqrt(1 + alpha**2 - 2*alpha*mu)
			ks_mu = k_lins * k_diff_ang
			del2 = del2_interp(ks_mu)
			mu_integrand[imu] = (del2 * smooth
									* polyee) / k_diff_ang**7
		mu_integral = trapz_2d(mus, mu_integrand)
		alpha_integrand[ialpha] = (del2_interp(alpha * k_lins) * smooth / alpha) * mu_integral
	alpha_integral = trapz_2d(alphas, alpha_integrand)
	P_ee_0 = pref * 2*pi**2 * alpha_integral / k_lins**3

	P_ee = grow_Pk(P_ee_0, z_lin, Dz, n=4)

	return P_ee, np.zeros_like(P_ee)
"""

def intrinsic_spectra(z_lin, k_lins, P_lin, z_growth, Dz):

	p_lin_z0=P_lin[0]
	#We need the growth function to rescale the power spectrum at the end...get growth at same z-values as P_lin
	assert np.isclose(z_growth[0],0.)
	assert ((z_growth[-1]>z_lin[-1]) or np.isclose(z_growth[-1],z_lin[-1]))
	if not np.allclose(z_lin, z_growth):
		Dz=interp.interp1d(z_growth,Dz)(z_lin)

	#Get del2(k) interpolator
	del2_interp = del4.Del2kInterpFunc(k_lins, p_lin_z0)
	#Compute prefactor from Mackey et al. 2001
	R=1
	sig2R=del2_interp.sigma2_r(R)
	pref = 5 * 0.4**2  / 32 / (sig2R**2)
	smooth=1.

	#Do the integral over alpha and mu
	alpha_integral = del4.alpha_mu_integral_mackey(k_lins, del2_interp, del4.g_ee)
	P_ee_0 = pref * 2*pi**2 * alpha_integral / k_lins**3
	#Get P_II(k,z) from P_II(k,0) using growth factor (It's D(z)^4, since we've integrated over P(k) twice)
	P_ee = del4.grow_Pk(P_ee_0, z_lin, Dz, n=4)

	#No GI, so return zeros for second term
	return P_ee, np.zeros_like(P_ee)

