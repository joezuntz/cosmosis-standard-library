''' 
	FASTPT is a numerical algorithm to calculate 
	1-loop contributions to the matter power spectrum or other 
	similar type integrals. 
	The method is presented in the paper 2016arXiv160304826M.
	Please cite this paper if you are using FASTPT
		
	J. E. McEwen (c) 2016 
	mcewen.24@osu.edu 
	
	The FASTPT class is the workhose of the FASTPT algorithm. 
	This class calculates integrals of the form 
	\int \frac{d^3q}{(2 \pi)^3} K(q,k-q) P(q) P(|k-q|) 
'''
from __future__ import division 
from __future__ import print_function


import numpy as np
from numpy.fft import fft, ifft , rfft, irfft , fftfreq
from numpy import exp, log, log10, cos, sin, pi, cosh, sinh , sqrt
from scipy.special import gamma 
import sys
from time import time 
from fastpt_extr import p_window, c_window, pad_left, pad_right
from matter_power_spt import P_13_reg 
from scipy.signal import convolve , fftconvolve 
from scipy.interpolate import interp1d
from gamma_funcs import g_m_vals, gamsn
from P_extend import k_extend 
from scipy.integrate import quad

log2=log(2.)
class FASTPT:
	
	def __init__(self,k,nu,param_mat=None,low_extrap=None,high_extrap=None,n_pad=None,verbose=False):

		
		# check for log spacing
		dk=np.diff(np.log(k))
		dk_test=np.ones_like(dk)*dk[0]
		
		log_sample_test='ERROR! FASTPT will not work if your k vector is not sampled evenly in log space!'
		
		np.testing.assert_array_almost_equal(dk, dk_test, decimal=4, err_msg=log_sample_test, verbose=False)
		
		# size of input array must be an even number 
		if (k.size % 2 != 0):
			raise ValueError('Input array must contain an even number of elements.')
			
		self.extrap=False		
		if (low_extrap is not None or high_extrap is not None):
			self.EK=k_extend(k,low_extrap,high_extrap)
			k=self.EK.extrap_k()
			self.extrap=True
			
		self.low_extrap=low_extrap
		self.high_extrap=high_extrap
		
		
		self.k_old=k
		
		delta_L=(log(np.max(k))-log(np.min(k)))/(k.size-1)
		
		if(n_pad !=None):
			self.id_pad=np.arange(k.size)+n_pad
			d_logk=delta_L
			k_pad=np.log(k[0])-np.arange(1,n_pad+1)*d_logk
			k_pad=np.exp(k_pad)
			k_left=k_pad[::-1]
			
			k_pad=np.log(k[-1])+np.arange(1,n_pad+1)*d_logk
			k_right=np.exp(k_pad)
			k=np.hstack((k_left,k,k_right))
			
			# check to make sure that the n padding sufficient to keep the 
			# FASTPT k_min less than 1/2 of the input k_min (added a plus 1 to be safe)
			n_pad_check=int(np.log(2)/delta_L) +1
			if (n_pad < n_pad_check): 
				print('Warning, you should consider increasing your zero padding to at least ', n_pad_check, ' .')
				print('So, that you ensure that k > 2k_min.')
				print(' k min in the FASTPT universe is ', k[0], ' while k min input is ', self.k_old[0])
						
	
		if(n_pad == None): 
			print('Your results are only good for k > 2k_min')
		
			
		# default parameters for standard P_22_reg
		if param_mat is None:
			param_mat=np.array([[0,0,0,0],[0,0,2,0],[0,0,4,0],[2,-2,2,0],\
							[1,-1,1,0],[1,-1,3,0],[2,-2,0,1] ])
	
		self.k=k
		self.k_size=k.size
		self.nu=nu
		self.p_mat=param_mat
		self.p_size=param_mat.shape[0]
		self.verbose=verbose
		self.n_pad=n_pad
		
		alpha=self.p_mat[:,0]
		beta=self.p_mat[:,1]
		l_Bessel=self.p_mat[:,2]
		type=self.p_mat[:,3]
		

		self.N=k.size
		
		# define eta_m and eta_n=eta_m
		omega=2*pi/(float(self.N)*delta_L)
		self.m=np.arange(-self.N//2,self.N//2+1) 
		self.eta_m=omega*self.m
		
		
		# define l and tau_l
		self.n_l=self.m.size + self.m.size - 1
		self.l=np.arange(-self.n_l//2+1,self.n_l//2+1)
		self.tau_l=omega*self.l
		
		#Q_m=np.zeros((param_mat.shape[0],self.N+1), dtype=complex) 
		self.pf=np.zeros((param_mat.shape[0])) 
		self.two_part_l=np.zeros((param_mat.shape[0],self.l.size), dtype=complex) 
		#Q_n=Q_m
		self.g_m=np.zeros((param_mat.shape[0],self.N+1), dtype=complex) 
		self.g_n=np.zeros((param_mat.shape[0],self.N+1), dtype=complex) 
		self.h_l=np.zeros((param_mat.shape[0],self.l.size),dtype=complex)
		
		self.p=-5-2*self.nu-alpha-beta
	
		for i in range(param_mat.shape[0]):
	
			sigma=l_Bessel[i]+1/2.
		
			# Define Q_m and Q_n and p 
			# use eta_m for Q_n, the value is the same 
			Q_m=3/2.+ nu + alpha[i] + 1j*self.eta_m
			Q_n=3/2.+ nu + beta[i] + 1j*self.eta_m
			p=-5-2*nu-alpha[i]-beta[i]
	
			self.g_m[i,:]=g_m_vals(sigma,Q_m)
			
			if (type[i]==1):
			
				# this is the special case, Corresponding to the regularized version of 
				# J_{2,-2,0,reg}(k)
				# get values for g_n 
				# again use eta_m. 
				s=2+nu + beta[i] 
				Q_n=s+ 1j*self.eta_m
		
				self.g_n[i,:]=gamsn(Q_n)
		
				#two_part_m=2**Q_m
				self.g_m[i,:]=self.g_m[i,:]*2.**Q_m
	
				# prefactor 
				self.pf[i]=(-1)**l_Bessel[i]/pi**3*np.sqrt(pi/2.)
			
				self.two_part_l[i,:]=np.ones(self.l.size)
			
			else:
				self.g_n[i,:]=g_m_vals(sigma,Q_n)
				# pre factor 
				self.pf[i]=(-1)**l_Bessel[i]/pi**2*2.**(2+2*nu+alpha[i]+beta[i])
			
				self.two_part_l[i,:]=exp(1j*self.tau_l*log2)
			
			# calculate h_l     
			#arg=(p+1-1j*tau_l)
			self.h_l[i,:]=gamsn(self.p[i]+1-1j*self.tau_l)
	
	def J_k(self,P,P_window=None,C_window=None):

		if(self.low_extrap is not None):
			P=self.EK.extrap_P_low(P)

		if(self.high_extrap is not None):
			P=self.EK.extrap_P_high(P)
			
		
		P_b=P*self.k_old**(-self.nu)
		
		if P_window is not None:
		# window the input power spectrum, so that at high and low k
		# the signal smoothly tapers to zero. This make the input
		# more like a periodic signal 
			
			if (self.verbose):
				print('windowing biased power spectrum')
			W=p_window(self.k_old,P_window[0],P_window[1])
			P_b=P_b*W 
			
		if (self.n_pad !=0): 
			P_b=np.pad(P_b, pad_width=(self.n_pad,self.n_pad), mode='constant', constant_values=0)
	
		c_m_positive=rfft(P_b)
		# End point should be divided by two. However, we always filter the Fourier coefficients (the last element is set to 
		# zero), so this really has no effect. 
		c_m_positive[-1]=c_m_positive[-1]/2.
		c_m_negative=np.conjugate(c_m_positive[1:])
		
		c_m=np.hstack((c_m_negative[::-1], c_m_positive))/float(self.N)
		
		if (C_window != None):
			# window the Fourier coefficients. 
			# This will damping the highest frequencies 
			
			if (self.verbose):
				print('windowing the Fourier coefficients')
			c_m=c_m*c_window(self.m,int(C_window*self.N//2.)) 
			
		
		A_out=np.zeros((self.p_size,self.k_size))
		for i in range(self.p_size):
	
			
			# convolve f_c and g_c 
			#C_l=np.convolve(c_m*self.g_m[i,:],c_m*self.g_n[i,:])
			C_l=fftconvolve(c_m*self.g_m[i,:],c_m*self.g_n[i,:])
	
			# multiply all l terms together 
			C_l=C_l*self.h_l[i,:]*self.two_part_l[i]
		
			# set up to feed ifft an array ordered with l=0,1,...,-1,...,N/2-1
			c_plus=C_l[self.l>=0]
			c_minus=C_l[self.l< 0]
		
			C_l=np.hstack((c_plus[:-1],c_minus))
			A_k=ifft(C_l)*C_l.size # multiply by size to get rid of the normalization in ifft
					
			A_out[i,:]=np.real(A_k[::2])*self.pf[i]*self.k**(-self.p[i]-2) 
			# note that you have to take every other element 
			# in A_k, due to the extended array created from the
			# discrete convolution 
		
		P_out=irfft(c_m[self.m>=0])*self.k**self.nu*float(self.N)
		if (self.n_pad !=0):
			# get rid of the elements created from padding 
			P_out=P_out[self.id_pad]
			A_out=A_out[:,self.id_pad] 
		
		return P_out, A_out
		
	def P22(self,P,P_window=None,C_window=None):
		
		Power, mat=self.J_k(P,P_window=P_window,C_window=C_window)  
		A=1219/1470.*mat[0,:]
		B=671/1029.*mat[1,:]
		C=32/1715.*mat[2,:]
		D=1/3.*mat[3,:]
		E=62/35.*mat[4,:]
		F=8/35.*mat[5,:]
		reg=1/3.*mat[6,:]

		return Power, 2*(A+B+C+D+E+F) + reg
		
	def one_loop(self,P,P_window=None,C_window=None):
	    
	    Ps,P22=self.P22(P,P_window,C_window)
	    P13=P_13_reg(self.k_old,Ps)
	    if (self.extrap):
	        _,P=self.EK.PK_original(P22+P13)
	        return P
	    
	    return P22+P13

	def P_bias(self,P,P_window=None,C_window=None): 
		# Quadraric bias Legendre components
		# See eg section B of Baldauf+ 2012 (arxiv: 1201.4827)
		# Note pre-factor convention is not standardized
		# Returns relevant correlations (including Wick contraction factors),
		# but WITHOUT bias values and other pre-factors.
		# Uses standard "full initialization" of J terms

		Power, mat=self.J_k(P,P_window=P_window,C_window=C_window)
		sig4=np.trapz(self.k_old**2*Power**2,x=self.k_old)/(2.*pi**2)
		#sig2=np.trapz(self.k_old**2*Power,x=self.k_old)/(2.*pi**2)
		
		Pd1d2=2.*(17./21*mat[0,:]+mat[4,:]+4./21*mat[1,:])
		Pd2d2=2.*(mat[0,:])
		Pd1s2=2.*(8./315*mat[0,:]+4./15*mat[4,:]+254./441*mat[1,:]+2./5*mat[5,:]+16./245*mat[2,:])
		Pd2s2=2.*(2./3*mat[1,:])
		Ps2s2=2.*(4./45*mat[0,:]+8./63*mat[1,:]+8./35*mat[2,:])

		if (self.extrap):
			_, Power=self.EK.PK_original(Power)
			_, Pd1d2=self.EK.PK_original(Pd1d2)
			_, Pd2d2=self.EK.PK_original(Pd2d2)
			_, Pd1s2=self.EK.PK_original(Pd1s2)
			_, Pd2s2=self.EK.PK_original(Pd2s2)
			_, Ps2s2=self.EK.PK_original(Ps2s2)

		return Power, Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig4
		

if __name__ == "__main__":
	# An example script to run FASTPT and get the plot for
	# P_22 + P_13 
	
	# load the data file 
	#d=np.loadtxt('Pk_Planck15.dat')
	d=np.loadtxt('Pk_test.dat') 
	# declare k and the power spectrum 
	k=d[:,0]; P=d[:,1]
	
	# set the parameters for the power spectrum window and
	# Fourier coefficient window 
	#P_window=np.array([.2,.2])  
	C_window=.75    
	
	# bias parameter and padding length 
	nu=-2; n_pad=1000
	
	from time import time
		
	# initialize the FASTPT class 
	# including extrapolation to higher and lower k  
	
	fastpt=FASTPT(k,nu,low_extrap=-5,high_extrap=3,n_pad=n_pad) 
	#fastpt=FASTPT(k,nu,n_pad=n_pad) 
	
	
	t1=time()   
	# without P_windowing (better if you are using zero padding) 
	P_spt=fastpt.one_loop(P,C_window=C_window) 
	#Ps,_,_,_,_,_,_=fastpt.P_bias(P,C_window=C_window) 
	
	t2=time()
	print('time', "%10.3f" %(t2-t1),'s')
	
	
	# make a plot 
	import matplotlib.pyplot as plt
	
	ax=plt.subplot(111)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylabel(r'$P(k)$', size=30)
	ax.set_xlabel(r'$k$', size=30)
	
	ax.plot(k,P,label='linear')
	ax.plot(k,P_spt, label=r'$P_{22}(k) + P_{13}(k)$' )
		
	plt.legend(loc=3) 
	plt.grid()
	plt.show()
