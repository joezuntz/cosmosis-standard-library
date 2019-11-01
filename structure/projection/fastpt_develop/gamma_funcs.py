''' This is the file that we keep all our Gamma function routines in.
	J.E. McEwen 
'''
import numpy as np
from numpy import exp, pi, sin, cos, log, sqrt 
from scipy.special import gamma 
	
def log_gamma(z):
	
	z=gamma(z)
	w=log(z)
	x=np.real(w)
	y=np.imag(w)
	return x,y
				

def g_m_vals(mu,q):

	imag_q= np.imag(q)
	
	g_m=np.zeros(q.size, dtype=complex)

	cut =200
	asym_q=q[np.absolute(imag_q) >cut]
	asym_plus=(mu+1+asym_q)/2.
	asym_minus=(mu+1-asym_q)/2.
	
	q_good=q[ (np.absolute(imag_q) <=cut) & (q!=mu + 1 + 0.0j)]

	alpha_plus=(mu+1+q_good)/2.
	alpha_minus=(mu+1-q_good)/2.
	
	g_m[(np.absolute(imag_q) <=cut) & (q!= mu + 1 + 0.0j)] =gamma(alpha_plus)/gamma(alpha_minus)

	#g_m[np.absolute(imag_q)>cut] = exp( (asym_plus-0.5)*log(asym_plus) - (asym_minus-0.5)*log(asym_minus) - asym_q )
	
	#g_m[np.absolute(imag_q)>cut] = exp( (asym_plus-0.5)*log(asym_plus) - (asym_minus-0.5)*log(asym_minus) - asym_q \
	#								+1./12 *(1./asym_plus - 1./asym_minus) +1./360.*(1./asym_minus**3 - 1./asym_plus**3) )
	
	# to higher order 								
	g_m[np.absolute(imag_q)>cut] = exp( (asym_plus-0.5)*log(asym_plus) - (asym_minus-0.5)*log(asym_minus) - asym_q \
	    +1./12 *(1./asym_plus - 1./asym_minus) +1./360.*(1./asym_minus**3 - 1./asym_plus**3) +1./1260*(1./asym_plus**5 - 1./asym_minus**5) )

	g_m[np.where(q==mu+1+0.0j)[0]] = 0.+0.0j
	
	return g_m
	
def gamsn(z):
	z=np.asarray(z, dtype=complex)
	result=sqrt(pi) /2. * 2**z *g_m_vals(0.5, z-0.5)
	return result
