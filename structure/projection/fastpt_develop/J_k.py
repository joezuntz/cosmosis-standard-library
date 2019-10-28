''' This file contains the routine to calculate 
	J_{\alpha, \beta, l}(k), as appears in 2.21 of the paper.
	
	It is the orginal FAST-PT code and has now been replaced by 
	FASTPT.py. 
	
	J. E. McEwen (c) 2016
''' 
from __future__ import division 
import numpy as np
from numpy.fft import fft, ifft , rfft, irfft , fftfreq
from numpy import exp, log, log10, cos, sin, pi, cosh, sinh , sqrt
from scipy.special import gamma 
import sys
from time import time 
from fastpt_extr import p_window, c_window, pad_left, pad_right
from gamma_funcs import g_m_vals, gamsn


log2=log(2.)
	

def check_conjugate(A):
	'''this function was only used for debugging. 
	It is not important for actually running FAST-PT
	This module is used to check that signal is real, 
	i.e. A(-omega)=A(omega)^*
	must be odd dimensional input 
	'''
	
	n=A.size
	
	for i in range(n/2):
		print('frequency indices'),  -n/2 + 1+ i, n/2-i
		print('values'), A[i], A[-i-1]
		print('sum of imaginary part'), np.imag(A[i] + A[-i-1])
		x=A[i] + A[-i-1]
		print('should be zero'), np.imag(x) 
	print('the zero frequency at '), n/2, '=', A[n/2], 
	return 1 
	
def J_k(k,P,param_matrix,nu=-2,P2=None,P_window=None, C_window=None,n_pad=500,verbose=False):

	# size of input array must be an even number 
	if (k.size % 2 != 0):
		raise ValueError('Input array must contain an even number of elements.')
	
	alpha=param_matrix[:,0]
	beta=param_matrix[:,1]
	l_Bessel=param_matrix[:,2]
	type=param_matrix[:,3]
	
	N=k.size
	
	delta_L=(log(np.max(k))-log(np.min(k)))/(N-1)
	
	P_b=P*k**(-nu)
	
	if P_window is not None:
		# window the input power spectrum, so that at high and low k
		# the signal smoothly tappers to zero. This make the input
		# more "like" a periodic signal 
		
		if (verbose):
			print('smoothing biased power spectrum')
		W=p_window(k,P_window[0],P_window[1])
		P_b=P_b*W 
	
	if (n_pad !=None):	
		# pad the edges. This helps with edge effects in Fourier space
		
		if (verbose):	
			print('padding the input signal with'), n_pad, 'zeros.'
		id_pad=np.arange(k.size)	
		k,P_b=pad_left(k,P_b,n_pad)	
		_,P=pad_left(k,P,n_pad)	
		k,P_b=pad_right(k,P_b,n_pad)
		_,P=pad_right(k,P,n_pad)	
		N=k.size
		id_pad=id_pad+n_pad
		
	
	# I take the real Fourier transform and then take the conjugate to
	# obtain the negative frequencies.  The convolution latter in the code requires 
	# the negative frequencies. 
	c_m_positive=rfft(P_b)
	c_m_negative=np.conjugate(c_m_positive[1:])
	c_m=np.hstack((c_m_negative[::-1], c_m_positive))/float(N)
		
	# frequency integers 
	n_c=c_m_positive.size
	m=np.arange(-n_c+1,n_c)
	
			
	# define eta_m and eta_n=eta_m
	omega=2*pi/(float(N)*delta_L)
	eta_m=omega*m
		
	# define l and tau_l
	n_l=c_m.size + c_m.size - 1
	l=l=np.arange(-n_l//2+1,n_l//2+1)
	tau_l=omega*l
	
	
	if (C_window != None):
		# window the Fourier coefficients. 
		# This will damping the highest frequencies 
		if (verbose):
			print('smoothing the Fourier coefficients')
		c_m=c_m*c_window(m,int(C_window*N//2.)) 
		
	# matrix for output 
	A_out=np.zeros((param_matrix.shape[0],k.size))
	for i in range(param_matrix.shape[0]):
	
		sigma=l_Bessel[i]+1/2.
		
		# Define Q_m and Q_n and p 
		# use eta_m for Q_n, the value is the same 
		Q_m=3/2.+ nu + alpha[i] + 1j*eta_m
		Q_n=3/2.+ nu + beta[i] + 1j*eta_m
		p=-5-2*nu-alpha[i]-beta[i]
	
		g_m=g_m_vals(sigma,Q_m)
			
		if (type[i]==1):
			
			# this is the special case, Corresponding to the regularized version of 
			# J_{2,-2,0,reg}(k)
			# get values for g_n 
			# again use eta_m. 
			s=2+nu + beta[i] 
			Q_n=s+ 1j*eta_m
		
			g_n=gamsn(Q_n)
		
			#two_part_m=2**Q_m
			g_m=g_m*2.**Q_m
	
			# prefactor 
			pf=(-1)**l_Bessel[i]/pi**3*np.sqrt(pi/2.)
			
			two_part_l=1
			
		else:
			g_n=g_m_vals(sigma,Q_n)
			# pre factor 
			pf=(-1)**l_Bessel[i]/pi**2*2.**(2+2*nu+alpha[i]+beta[i])
			
			two_part_l=exp(1j*tau_l*log2)

				
		# convolve f_c and g_c 
		C_l=np.convolve(c_m*g_m,c_m*g_n)

		# calculate h_l 	
		#arg=(p+1-1j*tau_l)
		h_l=gamsn(p+1-1j*tau_l)
				
		# multiply all l terms together 
		C_l=C_l*h_l*two_part_l
		
		# set up to feed ifft an array ordered with l=0,1,...,-1,...,N/2-1
		c_plus=C_l[l>=0]
		c_minus=C_l[l< 0]
		
		C_l=np.hstack((c_plus[:-1],c_minus))
		A_k=ifft(C_l)*C_l.size # multiply by size to get rid of the normalization in ifft
		
		
		A_out[i,:]=np.real(A_k[::2])*pf*k**(-p-2) # note that you have to take every other element 
												  # in A_k, due to the extended array created from the
												  # discrete convolution 
	P_out=irfft(c_m[m>=0])*k**nu*float(N)
	if (n_pad !=0):
		# get rid of the elements created from padding 
		P_out=P_out[id_pad]
		A_out=A_out[:,id_pad]	
		
	return P_out, A_out
	
		