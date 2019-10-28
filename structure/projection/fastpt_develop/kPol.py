import numpy as np
from J_table import J_table 
import sys

# The ordering is \alpha, \beta, l_1, l_2, l, A coefficient.
l_mat1= np.array([[-2,-2,0,0,0,2.],\
          [-2,-2,2,0,0,6.],\
          [-2,-2,0,0,2,1.],\
          [-2,-2,2,2,0,6.],\
          [-2,-2,1,1,1,-9.]], dtype=float)
          
l_mat2 = np.array([[-2,-2,0,0,0,1.],\
          [-2,-2,2,0,0,-2.],\
          [-2,-2,1,1,1,9.],\
          [-2,-2,2,2,0,-8.]], dtype=float)
          
l_mat3= np.array([[-2,-2,0,0,0,1.],\
          [-2,-2,2,0,0,-2.],
          [-2,-2,2,2,0,1.]], dtype=float)
          
def kPol():
	table1=np.zeros(10,dtype=float)
	for i in range(l_mat1.shape[0]):
	   
	    x=J_table(l_mat1[i])
	    table1=np.row_stack((table1,x))
	
	table2=np.zeros(10,dtype=float)
	for i in range(l_mat2.shape[0]):
	    
	    x=J_table(l_mat2[i])
	    table2=np.row_stack((table2,x))
		
	table3=np.zeros(10,dtype=float)
	for i in range(l_mat3.shape[0]):
	    
	    x=J_table(l_mat3[i])
	    table3=np.row_stack((table3,x))
	return table1[1:,:], table2[1:,:], table3[1:,:]

