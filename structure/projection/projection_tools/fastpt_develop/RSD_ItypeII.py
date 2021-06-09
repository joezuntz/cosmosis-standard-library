from __future__ import division 
import numpy as np
import sys
from time import time 
from numpy import log, sqrt, exp, pi
from scipy.signal import fftconvolve as convolve 


def P_Ap1(k,P,f):
	
	N=k.size
	n= np.arange(-N+1,N )
	dL=log(k[1])-log(k[0])
	s=n*dL
	cut=3
	high_s=s[s > cut]
	low_s=s[s < -cut]
	mid_high_s=s[ (s <= cut) &  (s > 0)]
	mid_low_s=s[ (s >= -cut) &  (s < 0)]


	Z1=lambda r : (3.*f *(6.*r - 22.*r**3 - 22. *r**5 + 6.*r**7 + 3.*log(np.absolute(r-1.)/(r+1.))*(-1. + r**2)**4) \
		-4.*r**2 *(9.*log(np.absolute(r-1.)/(r+1.))*(-1. + r**2)**3 + 2.*r *(19. - 24.* r**2 + 9.* r**4)))/r**3
	Z1_high=lambda r : -224.+(384-768*f/5.)*r**2+(-1152./5+2304./35*f)*r**4+(1152-256*f)/35.*r**6 + (128+34*f)/35.*r**8
	Z1_low=lambda r: 32./5-768*f/5. + (-1152+2304*f)/35./r**2 +(-128-256*f)/35./r**4 +(-264+34*f)/35./r**6 +(96./7-186.*f/35.)/r**8

	f_mid_low=Z1(exp(-mid_low_s))*exp(-mid_low_s)
	f_mid_high=Z1(exp(-mid_high_s))*exp(-mid_high_s)
	f_high = Z1_high(exp(-high_s))*exp(-high_s)
	f_low = Z1_low(exp(-low_s))*exp(-low_s)
	
	f=np.hstack((f_low,f_mid_low,-32.-96.*f,f_mid_high,f_high))
	# print(f)

	g= convolve(P, f) * dL
	g_k=g[N-1:2*N-1]
	ap1= 1./672./pi**2*k**2 * P*g_k 
	return ap1

def P_Ap3(k,P,f):
	
	N=k.size
	n= np.arange(-N+1,N )
	dL=log(k[1])-log(k[0])
	s=n*dL
	cut=3
	high_s=s[s > cut]
	low_s=s[s < -cut]
	mid_high_s=s[ (s <= cut) &  (s > 0)]
	mid_low_s=s[ (s >= -cut) &  (s < 0)]


	Z3=lambda r : f *(9 *log(np.absolute(r-1.)/(r+1.)) *(-1. + r**2)**3 *(-1. - 7.* r**2 + f *(-1. + r**2)) + 2*r *(9. - 185 *r**2 + 159.* r**4 - 63.* r**6 +\
       f *(9. - 33.* r**2 - 33.* r**4 + 9.* r**6)))/(r**3)
	Z3_high=lambda r : -448.*f+ (3072.*f-768*f**2)/5.*r**2+(-13824.*f+2304.*f**2)/35.*r**4+(2048.*f-256.*f**2)/35.*r**6 + (58*5*f+34*f**2)/35.*r**8
	Z3_low=lambda r: -704*f/5.-768*f**2/5. + 2304.*f**2 /35./r**2 +(-512*f-256*f**2)/35./r**4 +(-496*f+34*f**2)/35./r**6 +(774.*f-186.*f**2)/35./r**8

	f_mid_low=Z3(exp(-mid_low_s))*exp(-mid_low_s)
	f_mid_high=Z3(exp(-mid_high_s))*exp(-mid_high_s)
	f_high = Z3_high(exp(-high_s))*exp(-high_s)
	f_low = Z3_low(exp(-low_s))*exp(-low_s)
	
	f=np.hstack((f_low,f_mid_low,-160.*f-96.*f**2,f_mid_high,f_high))
	# print(f)

	g= convolve(P, f) * dL
	g_k=g[N-1:2*N-1]
	ap3= 1./672./pi**2*k**2 * P*g_k 
	return ap3

def P_Ap5(k,P,f):
	
	N=k.size
	n= np.arange(-N+1,N )
	dL=log(k[1])-log(k[0])
	s=n*dL
	
	cut=3
	high_s=s[s > cut]
	low_s=s[s < -cut]
	mid_high_s=s[ (s <= cut) &  (s > 0)]
	mid_low_s=s[ (s >= -cut) &  (s < 0)]


	Z5=lambda r : -f**2 *(9*log(np.absolute(r-1.)/(r+1.)) *(-1. + r**2)**3 *(1. + 3.* r**2) + 2.* r *(-9. + 109. *r**2 - 63.* r**4 + 27.* r**6))/(r**3)
	Z5_high=lambda r : -224.*f**2+(1152*f**2/5.)*r**2 - 1152.*f**2/7.*r**4+(128.*f**2)/5.*r**6 + (162.*f**2)/35.*r**8
	Z5_low=lambda r: -736.*f**2/5. + (1152*f**2)/35./r**2 -384.*f**2/35./r**4 -46.*f**2/7./r**6 +42.*f/5./r**8

	f_mid_low=Z5(exp(-mid_low_s))*exp(-mid_low_s)
	f_mid_high=Z5(exp(-mid_high_s))*exp(-mid_high_s)
	f_high = Z5_high(exp(-high_s))*exp(-high_s)
	f_low = Z5_low(exp(-low_s))*exp(-low_s)
	
	f=np.hstack((f_low,f_mid_low,-128.*f**2,f_mid_high,f_high))
	# print(f)

	g= convolve(P, f) * dL
	g_k=g[N-1:2*N-1]
	ap5= 1./672./pi**2*k**2 * P*g_k 
	return ap5

