''' This is the file that I keep the extra functions for FASTPT. 
	It mostly contains window functions and padding routines. 
	J. E. McEwen
'''

import numpy as np
from numpy import pi, cos, sin, log 
from numpy.fft import rfft, irfft, fft, ifft 
import sys
from scipy import signal
from scipy.special import erf
from scipy.signal import butter, lfilter, filtfilt, lfilter_zi

	
def p_window(k,log_k_left,log_k_right):
	
	log_k=np.log10(k)
	
	max=np.max(log_k)
	min=np.min(log_k)
	
	# all the log k to the left 
	# and to the right
	log_k_left=min+log_k_left
	log_k_right=max-log_k_right
		
	left=log_k[log_k <= log_k_left]
	right=log_k[log_k >= log_k_right]
	x_right=(right- right[right.size-1])/(right[0]-max)
	x_left=(min-left)/(min-left[left.size-1])
	
	W=np.ones(k.size)
	W[log_k <= log_k_left] = (x_left - 1/(2*pi)*sin(2*pi*x_left))
	W[log_k  >= log_k_right] = (x_right-  1/(2*pi)*sin(2*pi*x_right))
	
	return W 
	
	
def c_window(n,n_cut):

	n_right = n[-1] - n_cut
	n_left = n[0]+ n_cut 

	n_r=n[ n[:]  > n_right ] 
	n_l=n[ n[:]  <  n_left ] 
	
	theta_right=(n[-1]-n_r)/float(n[-1]-n_right-1) 
	theta_left=(n_l - n[0])/float(n_left-n[0]-1) 

	W=np.ones(n.size)
	W[n[:] > n_right]= theta_right - 1/(2*pi)*sin(2*pi*theta_right)
	W[n[:] < n_left]= theta_left - 1/(2*pi)*sin(2*pi*theta_left)
	
	return W
				
def pad_left(k,P,n_pad):
	d_logk=np.log10(k[1])-np.log10(k[0])
	
	k_pad=np.log10(k[0])-np.arange(1,n_pad+1)*d_logk
	k_pad=10**k_pad
	k_pad=k_pad[::-1]
	P_pad=np.zeros(n_pad)
	
	return np.hstack((k_pad,k)), np.hstack((P_pad,P))
	
def pad_right(k,P,n_pad):
	d_logk=np.log10(k[1])-np.log10(k[0])
	
	k_pad=np.log10(k[-1])+np.arange(1,n_pad+1)*d_logk
	k_pad=10**k_pad
	P_pad=np.zeros(n_pad)
	
	return np.hstack((k,k_pad)), np.hstack((P,P_pad))
	

def n_eff(k,P):
	ln_p=np.log(P)
	ln_k=np.log(k)
	
	return k[:-1], np.diff(ln_p)/np.diff(ln_k)


