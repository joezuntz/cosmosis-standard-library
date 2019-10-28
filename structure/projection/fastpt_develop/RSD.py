from __future__ import division 
import numpy as np
from J_table import J_table 
from get_nu1_nu2 import nu1_nu2
from copy import copy 
import sys

f=1
# f is an artifact from earlier work. It appears in the one of the tables in the paper, 
# we set f equal to 1, since it made our coding easier. 
l_mat_A = np.array([
          [-1,0,0,1,0,[68./21,2./3*f,0,0, 26./9*f,2./3*f**2,0,0, 10./63*f**2]],\
          [-1,0,2,1,0,[0,-68./21*f,0,0,340./63*f,-52./21*f**2,0,0,260./63*f**2]],\
          [-1,0,1,0,1,[2,124./35*f,0,0,-92./105*f,108./35*f**2, 0,0 ,-254./105*f**2]],\
          [-1,0,1,2,1,[0,-2*f, 0,0,10./3*f,-2*f**2,0,0,10./3*f**2]],\
          [-1,0,0,1,2,[16./21,4./3*f,0,0,4./9*f,4./3*f**2,0,0,-52./63*f**2]],\
          [-1,0,2,1,2,[0,-16./21*f,0,0,80./63*f,-32./21*f**2, 0,0,160./63*f**2]],\
          [-1,0,1,0,3,[0,16./35*f,0,0,-16./35*f,32./35*f**2,0,0,-32./35*f**2]],\
          [-2,1,1,0,0,[0,2./3*f,0, 0,-2./3*f,2./3*f**2,0,0,-2./3*f**2]],\
          [-2,1,0,1,1,[2.,0,0,0, 8./3*f,0,0,0 ,2./3*f**2]],\
          [-2,1,2,1,1,[0,-2*f, 0,0,10./3*f,-2*f**2,0,0,10./3*f**2]],\
          [-2,1,1,0,2,[0,4./3*f,0,0,-4./3*f,4./3*f**2,0,0,-4./3*f**2]]],dtype=object)
          
          
l_mat_B = np.array([
          [-1,-1,1,1,0,[-1./2,-3./10*f,-1./20*f**2,3./2,0,3./20*f**2,0,3./2*f,-21./20*f**2,0,0,131/100.*f**2]],\
          [-1,-1,3,1,0,[0,3./10*f,1./10*f**2,0,-3*f,-6./5*f**2,0,7./2*f,-3./10*f**2,0,0,47./25*f**2]],\
          [-1,-1,3,3,0,[0,0,-f**2/20.,0,0, 21./20*f**2,0,0,-63./20*f**2,0,0,231./100*f**2]],\
          [-1,-1,0,0,1,[1/2.,f/2.,5*f**2/36,-1./2,0,1./12*f**2,0,-f/2.,-f**2/12.,0,0,-5/36*f**2]],\
          [-1,-1,2,0,1,[0,-f/2.,-5./18*f**2,0,3*f,4./3*f**2,0,-5./2*f,f**2/6.,0,0,-11./9*f**2]],\
          [-1,-1,2,2,1,[0,0,5./36*f**2,0,0,-17./12*f**2,0,0,53./12*f**2,0,0,-113/36*f**2]]],dtype=object)


def RSD_vals(l_mat,id):
    idrow=[0,1,2,3,4,id]
    
    table=np.zeros(10,dtype=float)
    for i in range(l_mat.shape[0]):
        x=J_table(l_mat[i,idrow])
        
        table=np.row_stack((table,x))
    return table[1:,:] 


def RSDA():
    l_mat=l_mat_A
    x=RSD_vals(l_mat,5)
    y=copy(x[:,8])
    A=np.zeros((y.size,9))
    for i in range(y.size):
        A[i,0],A[i,1],A[i,2],A[i,3],A[i,4],A[i,5],A[i,6],A[i,7],A[i,8]=y[i]
    x[:,8]=np.ones(x.shape[0])  
    return x,A
    
def RSDB():
    l_mat=l_mat_B  
    x=RSD_vals(l_mat,5)
    y=copy(x[:,8])
    A=np.zeros((y.size,12))
    for i in range(y.size):
        A[i,0],A[i,1],A[i,2],A[i,3],A[i,4],A[i,5],A[i,6],A[i,7],A[i,8],A[i,9],A[i,10],A[i,11]=y[i]
    x[:,8]=np.ones(x.shape[0]) 
    return x,A
    
    
	
