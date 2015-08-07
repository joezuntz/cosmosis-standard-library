"""
The nonlinear matter power spectrum from e.g. The Cosmic Emulator has a limited
k-range. Extrapolate at the lower end using 
Delta^2(k) \propto k^(3+n_s) => P(k) \propto k^n_s (Mead et al. 2015) 
Extrapolate at the 
upper end using fitting funtion in eq. 5 of Harnois-Deraps:
P(k) \propto k^(alpha(z)-3.0), alpha = alpha_0*(1+z)^0.1, alpha_0=0.92
"""
from cosmosis.datablock import names, option_section
import numpy as np
import scipy.interpolate as interp

def setup(options):
	kmin = options.get_double(option_section, "kmin",1.e-5)
	kmax = options.get_double(option_section, "kmax",100.)
	alpha_0 = options.get_double(option_section, "alpha_0",0.92)
	dlog10k = options.get_double(option_section, "dlog10k", 0.03)
	print 'Interpolating/extrapolating P_nl over %f < k < %f, dlog10k=%f'%(kmin,kmax,dlog10k)
	return kmin,kmax,alpha_0,dlog10k

def execute(block, config):
	# load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
	nl = names.matter_power_nl
	cosmo = names.cosmological_parameters
	kmin,kmax,alpha_0,dlog10k = config
	z_nl,k_orig,p_nl=block.get_grid(nl,"z","k_h","p_k")

	extrap_lo,extrap_hi=True,True
	if k_orig[0]<kmin:
		if k_orig[-1]>kmax:
			return 0
		else:
			kmin=k_orig[0]
			extrap_lo=False
			
	if k_orig[-1]>kmax:
		kmax=k_orig[-1]
		extrap_hi=False

	logk_orig=np.log10(k_orig)
	log_p_interp=interp.RectBivariateSpline(z_nl,np.log10(k_orig),np.log10(p_nl))
	n_s = block[cosmo, "n_s"]

	logk_out=np.arange(np.log10(kmin),np.log10(kmax),dlog10k)
	low_k_inds=np.where(logk_out<logk_orig[0])[0]
	high_k_inds=np.where(logk_out>logk_orig[-1])[0]

	#Convenient to work with regular output grid...
	LOGK_out,Z_out=np.meshgrid(logk_out,z_nl)
	P_out_orig=10**log_p_interp(z_nl,logk_out,grid=True)
	P_out=P_out_orig.copy()
	if extrap_lo:
	    #Calculate P for lower and upper k ranges
	    #lower:
	    P_out_lower=(10**LOGK_out)**n_s
	    logk_match_lower=low_k_inds[-1]+1
	    ratio_lower=P_out_orig[:,logk_match_lower]/P_out_lower[:,logk_match_lower]
	    _,ratio_lower_grid=np.meshgrid(np.ones_like(logk_out),ratio_lower)
	    P_out_lower*=ratio_lower_grid	
	    lowk_inds=(LOGK_out<logk_orig[0])
	    P_out[lowk_inds]=P_out_lower[lowk_inds]

	#upper:
	if extrap_hi:
	    ALPHA=alpha_0*(1+Z_out)**0.1
	    P_out_higher=(10**LOGK_out)**(ALPHA-3.)
	    #Match P_out_lower and P_out_higher to P_out_orig at appropriate k

	    logk_match_higher=high_k_inds[0]-1
	    ratio_higher=P_out_orig[:,logk_match_higher]/P_out_higher[:,logk_match_higher]
	    _,ratio_higher_grid=np.meshgrid(np.ones_like(logk_out),ratio_higher)
	    P_out_higher*=ratio_higher_grid

	    #Now set lower and upper ranges of P_out_orig appropriately
	    high_k_inds=(LOGK_out>logk_orig[-1])
	    P_out[high_k_inds]=P_out_higher[high_k_inds]

	#Replace k and p_k
	k_out=10**logk_out
	block.replace_grid(nl, "z", z_nl, "k_h",k_out,"P_k", P_out)
	#test plots
	#import pylab
	#pylab.loglog(k_orig,p_nl[0],'b')
	#pylab.loglog(k_orig,p_nl[100],'b')
	#pylab.loglog(k_out,P_out_orig[0],'g')
	#pylab.loglog(k_out,P_out_orig[100],'g')	
	#pylab.loglog(k_out,P_out_lower[0],'r')
	#pylab.loglog(k_out,P_out_lower[100],'r')
	#pylab.loglog(k_out,P_out_higher[0],'m')
	#pylab.loglog(k_out,P_out_higher[100],'m')
	#pylab.loglog(k_out,P_out[0],'k+')
	#pylab.loglog(k_out,P_out[100],'k+')	
	#pylab.show()	
	#print logk_match_lower,logk_match_higher
	return 0

def cleanup(config):
	pass
