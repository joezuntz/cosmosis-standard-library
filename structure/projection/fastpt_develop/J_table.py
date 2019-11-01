''' This module tranforms a 
    set of l_1, l_2, and l in a 
    set of J_1, J_2, and J_k. 
    
''' 
import numpy as np
from numpy import pi 
import sys
from Wigner_symbols import three_j, six_j

def J_range(l_a,l_b):
	# returns J values for l_a, l_b such that
	# | l_a - l_b | <= J <= l_a + l_b
	start=np.absolute(l_a-l_b)
	end= l_a+l_b 
	
	length=end-start
	if length==0:
		return np.array([end])		
	J = np.arange(length+1)
	J += start
	return J 

def coeff_B(l1,l2,l,J1,J2,Jk):

	if ((J1+l2+l)%2==0)&((l1+J2+l)%2==0)&((l1+l2+Jk)%2==0):
		PF=(-1)**(l+(J1+J2+Jk)/2)  * (2*J1+1)*(2*J2+1)*(2*Jk+1) / pi**3
		result=PF *three_j(np.array([J1,l2,l]),np.array([0,0,0]) ) \
			  	  *three_j(np.array([l1,J2,l]),np.array([0,0,0]) ) \
			  	  *three_j(np.array([l1,l2,Jk]),np.array([0,0,0]) ) \
			  	  *three_j(np.array([J1,J2,Jk]),np.array([0,0,0]) ) \
			  	  *six_j(np.array([J1,J2,Jk,l1,l2,l]) ) 
	else:
		result = 0
	return result 

def J_table(params):
	alpha, beta, l1,l2,l, A=params
	
	
	# creates a table of J_1, J_2, and J_k values
	# from a provided set of l_1, l_2, l values 
	# the out put is in the form of 
	# alpha, beta, l1, l2, l, J1, J2, Jk, A, B

	table=np.zeros(10,dtype=float)
	J1_range = J_range(l2, l)
	J2_range = J_range(l1, l)
	Jk_range = J_range(l1, l2)
	for J1 in J1_range:
		for J2 in J2_range:
			for Jk in Jk_range:
				B = coeff_B(l1,l2,l,J1,J2,Jk)
				if (B!=0):		
					#x=np.array([alpha,beta,l1,l2,l,J1,J2,Jk,A,B],dtype=float)
					x=np.array([alpha,beta,l1,l2,l,J1,J2,Jk,A,B],dtype=object)
					table=np.row_stack((table,x))

	# return the table, excluding the first row of all zeros
	return table[1:,:]

