''' 
    This module extends the power spectrum down to lower and higher k, 
    by calculating the power law index at both ends of the array. The extension is to 
    help with edge effects and should be removed when returned from FASTPT.  
    
    J.E. McEwen (c) 2016
    mcewen.24@osu.edu 
'''

import numpy as np
from numpy import log, exp, log10 
import sys


class k_extend: 

    def __init__(self,k,low=None,high=None):
                
        self.DL=log(k[1])-log(k[0]) 
        
        if low is not None:
            if (low > log10(k[0])):
                low=log10(k[0])
                print('Warning, you selected a extrap_low that is greater than k_min. Therefore no extrapolation will be done.')
                #raise ValueError('Error in P_extend.py. You can not request an extension to low k that is greater than your input k_min.')
        
            low=10**low
            low=log(low)
            N=np.absolute(int((log(k[0])-low)/self.DL))
           
            if (N % 2 != 0 ):
                N=N+1 
            s=log(k[0]) -( np.arange(0,N)+1)*self.DL 
            s=s[::-1]
            self.k_min=k[0]
            self.k_low=exp(s) 
           
            self.k=np.append(self.k_low,k)
            self.id_extrap=np.where(self.k >=self.k_min)[0] 
            k=self.k
            

        if high is not None:
            if (high < log10(k[-1])):
                high=log10(k[-1])
                print('Warning, you selected a extrap_high that is less than k_max. Therefore no extrapolation will be done.')
                #raise ValueError('Error in P_extend.py. You can not request an extension to high k that is less than your input k_max.')
            
            high=10**high
            high=log(high)
            N=np.absolute(int((log(k[-1])-high)/self.DL))
            
            if (N % 2 != 0 ):
                N=N+1 
            s=log(k[-1]) + (np.arange(0,N)+1)*self.DL 
            self.k_max=k[-1]
            self.k_high=exp(s)
            self.k=np.append(k,self.k_high)
            self.id_extrap=np.where(self.k <= self.k_max)[0] 
            

        if (high is not None) & (low is not None):
            self.id_extrap=np.where((self.k <= self.k_max) & (self.k >=self.k_min))[0]
            
            
    def extrap_k(self):
        return self.k 
        
    def extrap_P_low(self,P):
      
        ns=(log(P[1])-log(P[0]))/self.DL
        Amp=P[0]/self.k_min**ns
        P_low=self.k_low**ns*Amp
        return np.append(P_low,P) 

    def extrap_P_high(self,P):
       
        ns=(log(P[-1])-log(P[-2]))/self.DL
        Amp=P[-1]/self.k_max**ns
        P_high=self.k_high**ns*Amp
        return np.append(P,P_high) 
    
    def PK_original(self,P): 
        return self.k[self.id_extrap], P[self.id_extrap]
    

    