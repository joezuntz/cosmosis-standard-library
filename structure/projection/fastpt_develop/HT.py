''' python version of python FFTLOG.
	Joseph E. McEwen
	McEwen Laboratories, copyright 2015 
	email: jmcewen314@gmail.com
'''
from __future__ import division 
from __future__ import print_function

import numpy as np
from numpy.fft import fft, ifft , fftshift, ifftshift , rfft, irfft 
from numpy import exp, log, log10, cos, sin, pi
from scipy.special import gamma 
#import fftlog # JAB's fftlog hack 
from time import time 
from numpy import gradient as grad
import sys

log2=log(2)

def asym_raised(k,log_k_left,log_k_right):

	log_k_left=2
	log_k_right=2
	print('side values', log_k_left, log_k_right)
	
	log_k=np.log10(k)
	
	max=np.max(log_k)
	min=np.min(log_k)
	
	log_k_left=min+log_k_left
	log_k_right=max-log_k_right
	
	print('side values', log_k_left, log_k_right )
	
	left=log_k[log_k <= log_k_left]
	right=log_k[log_k >= log_k_right]
	theta_right=(right[0]-right)/(max-right[0])*pi
	theta_left=(left[left.size-1]-left)/(min-left[left.size-1])*pi
	
	W=np.ones(k.size)
	W[log_k <= log_k_left] = (1+ cos(theta_left))/2.
	W[log_k  >= log_k_right] = (1+cos(theta_right))/2.
	return W 
	
	

def log_gamma(z):
	
	z=gamma(z)
	w=log(z)
	x=np.real(w)
	y=np.imag(w)
	return x,y
	
def get_k0(N,mu,q,r0,L,k0):
	
	
	kr=float(k0*r0)
	delta_L=L/N
	
	x=q + 1j*pi/delta_L
	
	x_plus=(mu+1+x)/2.
	x_minus=(mu+1-x)/2.
		
	rp,phip=log_gamma(x_plus)
	rm,phim=log_gamma(x_minus)
	
	arg=log(2/kr)/delta_L + (phip - phim)/pi 
	iarg=np.rint(arg)
	if ( arg != iarg):
		kr=kr*exp((arg-iarg)*delta_L)
		#kr=kr*exp((arg+iarg)*delta_L)		# Hamilton sign 
	
	return kr 
	
def u_m_vals(m,mu,q,kr,L):

	x=q + 1j*2*pi*m/L
	
	alpha_plus=(mu+1+x)/2.
	alpha_minus=(mu+1-x)/2.
		
	rp, phip=log_gamma(alpha_plus) 
	rm, phim=log_gamma(alpha_minus) 
	
	log_r=q*log2 + rp - rm 
	phi=2*pi*m/L*log(2./kr) + phip - phim 
	
	real_part=exp(log_r)*cos(phi)
	imag_part=exp(log_r)*sin(phi) 
	
	u_m=real_part + 1j*imag_part 
	
	# adjust endpoint, the N/2=m.size point 
	u_m[m.size-1]=np.real(u_m[m.size-1])
	return u_m
	
def fft_log(k,f_k,q,mu):

	if ((q+mu) < -1) :
		print('Error in reality condition for Bessel function integration.')
		print(' q+mu is less than -1.')
		print('See Abramowitz and Stegun. Handbook of Mathematical Functions pg. 486')
		sys.exit()
	
	if ( q > 1/2.) :
		print('Error in reality condition for Bessel function integration.')
		print(' q is greater than 1/2')
		print('See Abramowitz and Stegun. Handbook of Mathematical Functions pg. 486')
		sys.exit()
		
				
	N=f_k.size
	delta_L=(log(np.max(k))-log(np.min(k)))/(N-1)
	delta_L10=(np.log10(np.max(k))-np.log10(np.min(k)))/(N-1)
	L=(log(np.max(k))-log(np.min(k)))
		
	# find a better way to check if it is evenly spaced in log 
	diff=np.diff(np.log(k))
	diff=np.diff(diff)
	if (np.sum(diff) >=1e-7): #orig was 1e-10. this step should be updated.
		print('You need to send in data that is sampled evenly in logspace')
		print('total is',np.sum(diff))
		print('Terminating code in fft_log')
		sys.exit()
		
	#k_n=k0*exp(h_n)
	#h_n=n*delta_L
	# k_n=exp(log(k0))*exp(n*delta_L) 
	
	#log_k0=(log(np.max(k))+log(np.min(k)))/2.
	log_k0=log(k[N//2])
	k0=exp(log_k0)
	
	
	# Fourier transform input data 
	# and get m values, shifted so the zero point is at the center
	
	c_m=rfft(f_k)
	m=np.fft.rfftfreq(N,d=1.)*float(N)
	# make r vector 
	#kr=get_k0(float(N),mu,q,1/k0,L,k0)
	kr=1
	#print('this is good kr', kr)
	r0=kr/k0
	log_r0=log(r0)
	
	m=np.fft.rfftfreq(N,d=1.)*float(N)
	m_r=np.arange(-N//2,N//2)
	m_shift=np.fft.fftshift(m_r)
	#print(m_shift)
	
	#s-array 
	s=delta_L*(-m_r)+log_r0		
	id=m_shift
	r=10**(s[id]/log(10))
	
	#m_shift=np.fft.fftshift(m)
	
	# get h array 	
	h=delta_L*m + log_k0
			

	u_m=u_m_vals(m,mu,q,kr,L)

	b=c_m*u_m
		
	A_m=irfft(b)
	
	A=A_m[id]
	
	# reverse the order 
	A=A[::-1]
	r=r[::-1]
	
	if (q!=0):
		A=A*(r)**(-float(q))
		
	return r, A 

##########################################################################################
# End of fftlog algorithm 


def k_to_r(k,f_k,alpha_k, beta_r, mu, pf,q=0):
	t1=time()
	# module to calculate Hankel Transform
	# \int_0^\infty dk r A(k) J_mu(kr), via fftlog algorithm
	# Common application is for power spectrum:
	# \xi(r)= \int dk k^2 /(2 \pi^2) \sin(kr)/kr P(k) 
	# in which case 
	# alpha_k=1.5
	# beta_r=-1.5
	# mu=.5 
	# pf=(2*np.pi)**(-1.5)
	
	f_k=k**alpha_k*f_k
	
	r, A=fft_log(k,f_k,q,mu)

	f_r=pf*A*r**beta_r 
	t2=time()
	#print('time to run ', t2-t1 )
	return r, f_r 
	
def r_to_k(r,f_r,alpha_k, beta_r, mu, pf,q=0):
	t1=time()
	# module to calculate Hankel Transform
	# \int_0^\infty dr k A(r) J_mu(kr), via fftlog algorithm
	# Common application is for correlation function:
	# P(k)= 2 pi \int dr r^2  \sin(kr)/kr xi(r) 
	# in which case 
	# alpha_k=-1.5
	# beta_r=1.5
	# mu=.5 
	# pf=2 pi *sqrt(pi/2)
	
	f_r=r**beta_r*f_r
	#print('hello there')
	k, A=fft_log(r,f_r,q,mu)
	
	f_k=pf*A*k**alpha_k 
	t2=time()
	#print('time to run ', t2-t1)
	return k, f_k
	


if __name__=='__main__':
	# load Jonathan's power spectrum 	
	#data=np.loadtxt('Pk_fine.dat')
	data=np.loadtxt('Pk_Planck15.dat')
	#data=np.loadtxt('TSPT_out_z_1.5.dat')
	id=np.arange(data.shape[0])
	#id=id[::20]  # downsample factor 
	k=data[id,0]
	P=data[id,1]

	r,xi=k_to_r(k,P,1.5,-1.5,.5, (2*pi)**(-1.5))
	k2,ZZ=r_to_k(r,xi,-1.5,1.5,.5, 4*pi*np.sqrt(pi/2.))
	
	#X=fftlog.fftlog(P,k)
		
	import matplotlib.pyplot as plt

	ax=plt.subplot(221)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel(r'$r$', size=30)
	ax.set_ylabel(r'$\xi(r)$', size=30)
	#ax.set_xlim(60,200)

	# plot the correlation function 
	ax.plot(r,xi, color='black', label='fftlog python' )
#	ax.plot(X[:,0],X[:,1], '--', color='red', label='fftlog Fortran' )
	
	plt.legend()

	ax=plt.subplot(222)
	ax.set_xscale('log')
	ax.set_xlabel(r'$r$', size=30)
	ax.set_ylabel(r'ratio', size=30)
	ax.set_ylim(.9,1.1)
	

	# plot the ratio of the correlations between fortran and python methods
	#ax.plot(r,xi/X[:,1], color='black', label='fftlog python' )
	plt.grid()

	ax=plt.subplot(223)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel(r'$k$', size=30)
	ax.set_ylabel(r'$P(k)$', size=30)
	
	# plot original and roundtrip power spectrum 
	ax.plot(k2,ZZ, color='black', label='fftlog python' )
	ax.plot(k,P, '--', color='red',  )

	ax=plt.subplot(224)
	ax.set_xscale('log')
	ax.set_xlabel(r'$k$', size=30)
	ax.set_ylabel(r'log(ratio)', size=30)
	#ax.set_ylim(.99,1.01)

	#Delta=log(ZZ)-log(P)
	# plot the ratio of the original and roundtrip power spectrum 
	#ax.plot(k,Delta, color='red', label='fftlog python' )
	plt.grid()
	
	plt.show()
	
	#print(log(P)-log(P))
	#print(log(k2)-log(k))