from __future__ import division 
import numpy as np
from J_table import J_table 
import sys
from time import time 
from numpy import log, sqrt, exp, pi
from scipy.signal import fftconvolve as convolve 


def P_IA_B(k,P):
	
	N=k.size
	n= np.arange(-N+1,N )
	dL=log(k[1])-log(k[0])
	s=n*dL
	cut=3
	high_s=s[s > cut]
	low_s=s[s < -cut]
	mid_high_s=s[ (s <= cut) &  (s > 0)]
	mid_low_s=s[ (s >= -cut) &  (s < 0)]

	# For Zbar and take factor of 2 away
	Z1=lambda r : ((2.* r * (225.- 600.* r**2 + 1198.* r**4 - 600.* r**6 + 225.* r**8) + \
    				225.* (r**2 - 1.)**4 * (r**2 + 1.) * log(np.absolute(r-1)/(r+1)) )/(20160.* r**3) - 29./315*r**2 )/2.
	Z1_high=lambda r : ( (-16*r**4)/147. + (32*r**6)/441. - (37*r**8)/4704. - (19*r**10)/3528. + r**12/7056. + (5*r**14)/1176. - (5*r**16)/2016. )/2.
	Z1_low=lambda r: (-16./147 - 5/(2016.*r**12) + 5/(1176.*r**10) + 1/(7056.*r**8) - 19/(3528.*r**6) - 37/(4704.*r**4) + 32/(441.*r**2) )/2.

	# # 
	# Z1 = Z1/2.
	# Z1_high = Z1_high/2.
	# Z1_low = Z1_low/2.

	f_mid_low=Z1(exp(-mid_low_s))*exp(-mid_low_s)
	f_mid_high=Z1(exp(-mid_high_s))*exp(-mid_high_s)
	f_high = Z1_high(exp(-high_s))*exp(-high_s)
	f_low = Z1_low(exp(-low_s))*exp(-low_s)
	
	f=np.hstack((f_low,f_mid_low,-1./42.,f_mid_high,f_high))
	# print(f)

	# For Z without subtracting low-k limit
	# Z1=lambda r : (2.* r * (225.- 600.* r**2 + 1198.* r**4 - 600.* r**6 + 225.* r**8) + \
 #    				225.* (r**2 - 1.)**4 * (r**2 + 1.) * log(np.absolute(r-1)/(r+1)) )/(20160.* r**3)
	# Z1_high=lambda r : (29*r**2)/315. - (16*r**4)/147. + (32*r**6)/441. - (37*r**8)/4704. - (19*r**10)/3528. + r**12/7056. + (5*r**14)/1176. - (5*r**16)/2016.
	# Z1_low=lambda r: -16./147- 5/(2016.*r**12) + 5/(1176.*r**10) + 1/(7056.*r**8) - 19/(3528.*r**6) - 37/(4704.*r**4) + 32/(441.*r**2) + (29*r**2)/315.

	# f_mid_low=Z1(exp(-mid_low_s))*exp(-mid_low_s)
	# f_mid_high=Z1(exp(-mid_high_s))*exp(-mid_high_s)
	# f_high = Z1_high(exp(-high_s))*exp(-high_s)
	# f_low = Z1_low(exp(-low_s))*exp(-low_s)
	
	# f=np.hstack((f_low,f_mid_low,2./45,f_mid_high,f_high))

	g= convolve(P, f) * dL
	g_k=g[N-1:2*N-1]
	IA_B= k**3/(2.*pi**2) * P*g_k 
	return IA_B


def IA_A():
	# Ordering is \alpha, \beta, l_1, l_2, l, A coeficient 
	l_mat_IAA= np.array([[0,0,0,0,0,-31./210],\
          [0,0,2,0,0,-34./63],\
          [0,0,0,0,2,-47./147],\
          [0,0,2,0,2,-8./63],\
          [0,0,1,1,1,93./70],\
          [0,0,1,1,3,6./35],\
          [0,0,0,0,4,-8./245],\
          [1,-1,0,0,1,-3./10],\
          [1,-1,2,0,1,-1./3],\
          [1,-1,1,1,0,1./2],\
          [1,-1,1,1,2,1.],\
          [1,-1,0,2,1,-1./3],\
          [1,-1,0,0,3,-1./5]], dtype=float)
	table=np.zeros(10,dtype=float)
	for i in range(l_mat_IAA.shape[0]):
		x=J_table(l_mat_IAA[i])
		table=np.row_stack((table,x))
	return table[1:,:] 
	
def IA_DEE():
	l_mat_DEE= np.array([[0,0,0,0,0,-43./540],\
          [0,0,2,0,0,-167./756],\
          [0,0,4,0,0,-19./105],\
          [0,0,0,0,2,1./18],\
          [0,0,2,0,2,-7./18],\
          [0,0,1,1,1,11./20],\
          [0,0,3,1,1,19./20],\
          [0,0,2,2,0,-19./54]], dtype=float)
	table=np.zeros(10,dtype=float)
	for i in range(l_mat_DEE.shape[0]):
		x=J_table(l_mat_DEE[i])
		table=np.row_stack((table,x))
	return table[1:,:] 

def IA_DBB():
	l_mat_DBB= np.array([[0,0,0,0,0,13./135],\
          [0,0,2,0,0,86./189],\
          [0,0,4,0,0,16./105],\
          [0,0,0,0,2,2./9],\
          [0,0,2,0,2,4./9],\
          [0,0,1,1,1,-13./15],\
          [0,0,3,1,1,-4./5],\
          [0,0,2,2,0,8./27]], dtype=float)
	table=np.zeros(10,dtype=float)
	for i in range(l_mat_DBB.shape[0]):
		x=J_table(l_mat_DBB[i])
		table=np.row_stack((table,x))
	return table[1:,:] 
