import numpy as np
from J_table import J_table 
import sys

# Ordering is \alpha, \beta, l_1, l_2, l, A coeficient 
l_mat_E= np.array([[0,0,0,0,0,16./81],\
		[0,0,2,0,0,713./1134],\
		[0,0,4,0,0,38./315],\
		[0,0,2,2,0,95./162],\
		[0,0,1,1,1,-107./60],\
		[0,0,3,1,1,-19./15],\
		[0,0,0,0,2,239./756],\
		[0,0,2,0,2,11./9],\
		[0,0,2,2,2,19./27],\
		[0,0,1,1,3,-7./10],\
		[0,0,0,0,4,3./35]],dtype=float)

l_mat_B= np.array([[0,0,0,0,0,-41./405],\
		[0,0,2,0,0,-298./567],\
		[0,0,4,0,0,-32./315],\
		[0,0,2,2,0,-40./81],\
		[0,0,1,1,1,59./45],\
		[0,0,3,1,1,16./15],\
		[0,0,0,0,2,-2./9],\
		[0,0,2,0,2,-20./27],\
		[0,0,2,2,2,-16./27],\
		[0,0,1,1,3,2./5]],dtype=float)
def IA_tt():
	tableE=np.zeros(10,dtype=float)
	for i in range(l_mat_E.shape[0]):
		xE=J_table(l_mat_E[i])
		tableE=np.row_stack((tableE,xE))

	tableB=np.zeros(10,dtype=float)
	for i in range(l_mat_B.shape[0]):
		xB=J_table(l_mat_B[i])
		tableB=np.row_stack((tableB,xB))
	
	return tableE[1:,:] , tableB[1:,:]
	
