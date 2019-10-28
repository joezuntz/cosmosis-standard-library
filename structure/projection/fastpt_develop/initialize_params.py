''' 
	This is the file that we keep all our functions that only 
	depend on the k-grid. These routines are used to set up the initialization 
	quatities for each FAST-PT parameter matrix.
'''
import numpy as np
from numpy import exp, pi, sin, cos, log, sqrt 
from scipy.special import gamma 
from get_nu1_nu2 import nu1_nu2
import sys 
	
log2=log(2.)	
	
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
	
	g_m[np.absolute(imag_q)>cut] = exp( (asym_plus-0.5)*log(asym_plus) - (asym_minus-0.5)*log(asym_minus) - asym_q \
									+1./12 *(1./asym_plus - 1./asym_minus) +1./360.*(1./asym_minus**3 - 1./asym_plus**3) )

	g_m[np.where(q==mu+1+0.0j)[0]] = 0.+0.0j
	
	return g_m
	
def gamsn(z):
	z=np.asarray(z, dtype=complex)
	result=sqrt(pi) /2. * 2.**z *g_m_vals(0.5, z-0.5)
	return result

def scalar_stuff(p_mat,nu,N,m,eta_m,l, tau_l):

	alpha=p_mat[:,0]
	beta=p_mat[:,1]
	
	l_Bessel=p_mat[:,2]
	type=p_mat[:,3]
	
	pf=np.zeros((p_mat.shape[0]))
	two_part_l=np.zeros((p_mat.shape[0],l.size), dtype=complex) 
	g_m=np.zeros((p_mat.shape[0],N+1), dtype=complex) 
	g_n=np.zeros((p_mat.shape[0],N+1), dtype=complex) 
	h_l=np.zeros((p_mat.shape[0],l.size),dtype=complex)
	p=np.zeros(p_mat.shape[0])
	for i in range(p_mat.shape[0]):
		
		sigma=l_Bessel[i]+1/2.
		
		# Define Q_m and Q_n and p 
		# use eta_m for Q_n, the value is the same 
		Q_m=3/2.+ nu + alpha[i] + 1j*eta_m
		Q_n=3/2.+ nu + beta[i] + 1j*eta_m
		p[i]=-5-2*nu-alpha[i]-beta[i]
		
		
		g_m[i,:]=g_m_vals(sigma,Q_m)
		
		if (type[i]==1):

			# this is the special case, Corresponding to the regularized version of 
			# J_{2,-2,0,reg}(k)
			# get values for g_n 
			# again use eta_m. 
			s=2+nu + beta[i] 
			Q_n=s+ 1j*eta_m
			
			g_n[i,:]=gamsn(Q_n)
			
			#two_part_m=2**Q_m
			g_m[i,:]=g_m[i,:]*2.**Q_m
			
			# prefactor
			pf[i]=(-1)**l_Bessel[i]/pi**3*np.sqrt(pi/2.)
			two_part_l[i,:]=np.ones(l.size)
	
		else:
			g_n[i,:]=g_m_vals(sigma,Q_n)
			
			# prefactor
	
			pf[i]=(-1)**l_Bessel[i]/pi**2*2.**(2+2*nu+alpha[i]+beta[i])
			two_part_l[i,:]=exp(1j*tau_l*log2)
			
		# calculate h_l
		h_l[i,:]=gamsn(p[i]+1-1j*tau_l)
		
	return pf,p, g_m, g_n, two_part_l, h_l 
		
def tensor_stuff(p_mat,N,m,eta_m,l, tau_l):

	alpha=p_mat[:,0]
	beta=p_mat[:,1]
	
	nu1,nu2=nu1_nu2(alpha,beta)
	
	J_1=p_mat[:,2]
	J_2=p_mat[:,3]
	J=p_mat[:,4]
	A=p_mat[:,5]
	B=p_mat[:,6]
	
	pf=np.zeros(p_mat.shape[0])      
	g_m=np.zeros((p_mat.shape[0],N+1), dtype=complex) 
	g_n=np.zeros((p_mat.shape[0],N+1), dtype=complex) 
	h_l=np.zeros((p_mat.shape[0],l.size),dtype=complex)
	
	p=3.+ nu1+ nu2 +alpha + beta

	for i in range(p_mat.shape[0]):
		
		sigma_1=J_1[i]+1/2.
		sigma_2=J_2[i]+1/2.          
		sigma_3=J[i]+1/2.
	
		# Define Q_m and Q_n and p 
		# use eta_m for Q_n, the value is the same 
		# Q_m=3/2.+ self.nu1[i] + alpha[i] + 1j*self.eta_m
		# Q_n=3/2.+ self.nu2[i] + beta[i] + 1j*self.eta_m
		
		Q_m=3/2.+ nu1[i] + alpha[i] + 1j*eta_m
		Q_n=3/2.+ nu2[i] + beta[i] + 1j*eta_m         

		#p=3. + nu1+nu2 + alpha[i] + beta[i]

		g_m[i,:]=g_m_vals(sigma_1,Q_m)
		g_n[i,:]=g_m_vals(sigma_2,Q_n)

		pf[i] = A[i] * B[i] * pi**(1.5) /8.
		
		
		# calculate h_l     
		#arg=(p+1-1j*tau_l)
		h_l[i,:]=g_m_vals(sigma_3, -1.5 - p[i] - 1j*tau_l)
	
	return pf, p, nu1, nu2, g_m, g_n, h_l