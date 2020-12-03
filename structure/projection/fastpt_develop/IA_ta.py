from __future__ import division 
import numpy as np
from J_table import J_table 
import sys
from time import time 
from numpy import log, sqrt, exp, pi
from scipy.signal import fftconvolve as convolve 


def P_IA_deltaE2(k,P):
	
	N=k.size
	n= np.arange(-N+1,N )
	dL=log(k[1])-log(k[0])
	s=n*dL
	cut=3
	high_s=s[s > cut]
	low_s=s[s < -cut]
	mid_high_s=s[ (s <= cut) &  (s > 0)]
	mid_low_s=s[ (s >= -cut) &  (s < 0)]

	# For Zbar
	Z1=lambda r : 30. + 146*r**2 - 110*r**4 + 30*r**6 + log(np.absolute(r-1.)/(r+1.))*(15./r - 60.*r + 90*r**3 - 60*r**5 + 15*r**7)
	Z1_high=lambda r : 256*r**2 - 256*r**4 + (768*r**6)/7. - (256*r**8)/21. + (34*r**10)/21. - (62*r**12)/7. + (190*r**14)/21. - (10*r**16)/3.
	Z1_low=lambda r: 768./7 - 10/(3.*r**10) + 190/(21.*r**8) - 62/(7.*r**6) + 34/(21.*r**4) - 256/(21.*r**2)

	f_mid_low=Z1(exp(-mid_low_s))*exp(-mid_low_s)
	f_mid_high=Z1(exp(-mid_high_s))*exp(-mid_high_s)
	f_high = Z1_high(exp(-high_s))*exp(-high_s)
	f_low = Z1_low(exp(-low_s))*exp(-low_s)
	
	f=np.hstack((f_low,f_mid_low,96.,f_mid_high,f_high))
	# print(f)

	g= convolve(P, f) * dL
	g_k=g[N-1:2*N-1]
	deltaE2= k**3/(896.*pi**2) * P*g_k 
	return deltaE2


def IA_deltaE1():
	# Ordering is \alpha, \beta, l_1, l_2, l, A coeficient 
	l_mat_deltaE1= np.array([[0,0,0,2,0,17./21],\
          [0,0,0,2,2,4./21],\
          [1,-1,0,2,1,1./2],\
          [-1,1,0,2,1,1./2]], dtype=float)
	table=np.zeros(10,dtype=float)
	for i in range(l_mat_deltaE1.shape[0]):
		x=J_table(l_mat_deltaE1[i])
		table=np.row_stack((table,x))
	return table[1:,:] 



def IA_0E0E():
	# Ordering is \alpha, \beta, l_1, l_2, l, A coeficient 
	l_mat_0E0E= np.array([[0,0,0,0,0,29./90],\
          [0,0,2,0,0,5./63],\
          [0,0,2,2,0,19./18],\
          [0,0,0,4,0,19./35]], dtype=float)
	table=np.zeros(10,dtype=float)
	for i in range(l_mat_0E0E.shape[0]):
		x=J_table(l_mat_0E0E[i])
		table=np.row_stack((table,x))
	return table[1:,:] 


def IA_0B0B():
	# Ordering is \alpha, \beta, l_1, l_2, l, A coeficient 
	l_mat_0B0B= np.array([[0,0,0,0,0,2./45],\
          [0,0,2,0,0,-44./63],\
          [0,0,2,2,0,-8./9],\
          [0,0,0,4,0,-16./35],\
          [0,0,1,1,1,2.]], dtype=float)
	table=np.zeros(10,dtype=float)
	for i in range(l_mat_0B0B.shape[0]):
		x=J_table(l_mat_0B0B[i])
		table=np.row_stack((table,x))
	return table[1:,:] 
