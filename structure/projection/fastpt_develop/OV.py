import numpy as np
from J_table import J_table 
import sys

# Ordering is \alpha, \beta, l_1, l_2, l, A coeficient 
l_mat= np.array([[-2,0,0,0,0,2./3],\
          [-2,0,0,2,0,-2./3],\
          [0,-2,0,0,0,-2./3],\
          [0,-2,0,2,0,2./3]], dtype=float)
def OV():
	table=np.zeros(10,dtype=float)
	for i in range(l_mat.shape[0]):
		x=J_table(l_mat[i])
		table=np.row_stack((table,x))
		
	return table[1:,:] 
	
