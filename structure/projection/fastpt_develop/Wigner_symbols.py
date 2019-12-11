''' Wigner symbols. 
	J.E. McEwen 2016
	
	Algorithm based on the sympy implimentation of sympy.physics.wigner, 
	which was based off of:
	
		[Rasch03] J. Rasch and A. C. H. Yu, 'Efficient Storage Scheme for
		Pre-calculated Wigner 3j, 6j and Gaunt Coefficients', SIAM
		J. Sci. Comput. Volume 25, Issue 4, pp. 1416-1428 (2003)
	
	Additional routines have been added for special configurations of Wigner 3 j symbols
	
	The code has not be written in a way to handle large input values; beware of precission issues
	associated with large values aquired from factorials. 
'''

import numpy as np
try:
	from scipy.misc import factorial
except ImportError:
	from scipy.special import factorial
import sys


# store factorials in an array
def factorial_list(N):
	x=np.arange(N+1)
	return factorial(x) 

# old way for the factorials, no longer used. 
# def factorial_list(N):
# 
# 	FL=[1]
# 	# N = highest factorial needed for computation 
# 	# note_to_self : maybe it is faster to use the scipy factorial 
# 	# function and then append to Factorial_list
# 	
# 	#if N > len(FL): 
# 	if N > 1: 
# 		for i in range(len(FL), int(N+1)):
# 			FL.append(FL[i-1]*i)
# 		return FL[:int(N) + 1]
# 	else: 
# 		return FL
	

def three_j(j,m): 

	j_1,j_2,j_3=j
	m_1,m_2,m_3=m 

	# symmetry conditions 
	if int(j_1 * 2) != j_1 * 2 or int(j_2 * 2) != j_2 * 2 or \
			int(j_3 * 2) != j_3 * 2:
		raise ValueError("j values must be integer or half integer, error in three_j)")
	if int(m_1 * 2) != m_1 * 2 or int(m_2 * 2) != m_2 * 2 or \
			int(m_3 * 2) != m_3 * 2:
		raise ValueError("m values must be integer or half integer, error in three_m")
	if m_1 + m_2 + m_3 != 0:
		return 0
	

	PF= np.int((-1) ** int(j_1 - j_2 - m_3))
	M=-m_3; 
	
	a1 = j_1 + j_2 - j_3
	if a1 < 0:
		return 0
	a2 = j_1 - j_2 + j_3
	if a2 < 0:
		return 0
	a3 = -j_1 + j_2 + j_3
	if a3 < 0:
		return 0
	if (abs(m_1) > j_1) or (abs(m_2) > j_2) or (abs(m_3) > j_3):
		return 0
		
	# special case identities, taken from the Mathematica website for 3-j symbols 
	if ( (j_1==j_2) & (j_3==0) & (m_1==-m_2) & (m_3==0) ):
		
		return (-1)**(j_1-m_1)/np.sqrt(2*j_1+1) 
		
	#if ( (m_1==0) & (m_2==0) & (m_3==0) 	
	
	
	max_factor=max(j_1 + j_2 + j_3 + 1, j_1 + abs(m_1), j_2 + abs(m_2),
				  j_3 + abs(m_3))
	
	FL=factorial_list(max_factor) 
	
	Sqrt_Arg=(FL[int(j_1+j_2-j_3)] \
			 *FL[int(j_1-j_2+j_3)] \
			 *FL[int(-j_1+j_2+j_3)] \
			 *FL[int(j_1-m_1)] \
			 *FL[int(j_1+m_1)] \
			 *FL[int(j_2-m_2)] \
			 *FL[int(j_2+m_2)] \
			 *FL[int(j_3-m_3)] \
			 *FL[int(j_3+m_3)] )/FL[int(j_1+j_2+j_3+1)]
			 
	Sqrt_part=np.sqrt(Sqrt_Arg)
	
	# need to fix this 
	#if Sqrt_part.is_complex:
	#	Sqrt_part=Sqrt_part.as_real_imag()[0]
	
	i_min = max(-j_3 + j_1 + m_2, -j_3 + j_2 - m_1, 0)
	i_max = min(j_2 + m_2, j_1 - m_1, j_1 + j_2 - j_3)
	
	Sum=0
	for i in range(int(i_min),int(i_max) + 1): 
		denom=FL[i] \
			*FL[int(i+j_3-j_1-m_2)] \
			*FL[int(j_2 + m_2- i)] \
			*FL[int(j_1-i-m_1)] \
			*FL[int(i  + j_3-j_2 + m_1)] \
			*FL[int(j_1 +j_2 -j_3 -i)]
		Sum=Sum+np.int((-1)**i)/float(denom) 
	
	
	# might have to reset FL 
	return Sum*Sqrt_part*PF 

def Delta_coef(a,b,c,prec=None):
	
	if int(a + b- c) != (a + b - c):
		raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
	if int(a + c - b) != (a + c - b):
		raise ValueError("j values must be integer or half integer and fulfill the triangle relation")	
	if int(b + c - a) != (b + c - a):
		raise ValueError("j values must be integer or half integer and fulfill the triangle relation")
	if (a+ b - c) < 0:
		return 0
	if (a + c - b) < 0:
		return 0
	if (b + c - a) < 0:
		return 0

	max_factor=max(a + b -c, a+c-b, b+c -a, a+b +c + 1) 
	
	FL=factorial_list(max_factor) 
	
	Sqrt_Arg=float(FL[int(a + b -c)] \
			  *FL[int(a + c - b)] \
			  *FL[int(b + c - a)])/float(FL[int(a + b + c + 1)])
	
	Sqrt_part=np.sqrt(Sqrt_Arg)
	#if prec:
	#		Sqrt_part=Sqrt_part.evalf(
	
	return Sqrt_part
	 
def Racah(a,b,c,d,e,f, prec=None):
	# the Racah symbol 
		
	PF=Delta_coef(a,b,e,prec) \
		*Delta_coef(c,d,e,prec) \
		*Delta_coef(a,c,f,prec) \
		*Delta_coef(b,d,f, prec) 
	
	if PF==0:
		return 0 
	
	i_min = max(a + b + e, c + d + e, a + c + f, b + d + f)
	i_max = min(a + b+ c + d, a + d + e + f, b + c + e + f)
	
	max_factor=max(i_max + 1, a + b + c + d, a + d + e + f,b + c + e + f)
	
	FL=factorial_list(max_factor) 
	
	Sum=0
	for i in range(int(i_min), int(i_max)+ 1): 
		
		denom=FL[int(i-a-b-e)]\
			*FL[int(i-c-d-e)]\
			*FL[int(i-a-c-f)]\
			*FL[int(i-b-d-f)]\
			*FL[int(a + b + c + d - i)]\
			*FL[int(a + d + e + f - i)]\
			*FL[int(b + c + e + f - i)] 
		
		Sum=Sum+((-1)**i*FL[i+1])/float(denom)
	
	return PF*Sum*(-1)**int(a+b+c+d) 
	
def six_j(j): 
	j_1,j_2,j_3,j_4,j_5,j_6=j 
	
	return (-1)**int(j_1+j_2+j_4 +j_5)*Racah(j_1, j_2, j_5, j_4, j_3, j_6,)


if __name__=="__main__":	
	# j=np.array([2,6,4])
# 	m=np.array([0,0,0]) 
# 	j_1,j_2,j_3=j
# 	m_1,m_2,m_3=m
# 	W=wigner.wigner_3j(j_1,j_2,j_3, m_1,m_2,m_3)
# 	
# 	print(W
# 	
# 	print(three_j(j,m)
# 	j=np.array([3,3,3,3,3,3]) 
# 	print('my six j',six_j(j) 
# 	W=wigner.wigner_6j(3,3,3,3,3,3)
# 	print('6 j',W, -1/14.
# 
# 	print(six_j(np.array([5,5,5,5,5,5])), 1/52.

	print('some test cases for the 3-j symbols')
	print('test 1')
	print('----------------------------------------')
	print('j=1,2,3 & m=0,0,0 => ', -np.sqrt(3/35.))
	j=np.array([1,2,3]); m=np.array([0,0,0]) 
	print('check' , three_j(j,m))
	
	print('test 2')
	print('----------------------------------------')
	print('j=4,5,9 & m=0,0,0 => ', -21*np.sqrt(2/46189.))
	j=np.array([4,5,9]); m=np.array([0,0,0]) 
	print('check' , three_j(j,m))
	
	print('test 3')
	print('----------------------------------------')
	print('j=4,5,6 & m=1,0,-1 => ', -2*np.sqrt(1/429.))
	j=np.array([4,5,6]); m=np.array([1,0,-1]) 
	print('check' , three_j(j,m))
	
	print('some test cases for the 6-j symbols')
	print('test 1')
	print('----------------------------------------')
	print('j=4,4,4,4,4,4 => ', -467/18018.)
	j=np.array([4,4,4,4,4,4])
	print('check' , six_j(j))
	
	print('test 2')
	print('----------------------------------------')
	print('j=1,1,1,1,1,1 => ', 1/6.)
	j=np.array([1,1,1,1,1,1])
	print('check' , six_j(j))
	
	print('test 3')
	print('----------------------------------------')
	print('j=1,2,3,1,3,1 => ', 1/5.)
	j=np.array([1,2,3,1,2,1])
	print('check' , six_j(j))