''' This is the script to make figure 1. of the paper
	This script is a a good exampel of how to implement 1-loop calculations.
	See line 24 (or around line 24 ) for the call to FAST-PT
	J. E. McEwen
	email: jmcewen314@gmail.com
	
'''

import numpy as np 
#from matter_power_spt import one_loop
import FASTPT 
from time import time 

# load the input power spectrum data 
d=np.loadtxt('Pk_test.dat')
k=d[:,0]
P=d[:,1]

# use if you want to interpolate data 
#from scipy.interpolate import interp1d 
#power=interp1d(k,P)
#k=np.logspace(np.log10(k[0]),np.log10(k[-1]),3000)
#P=power(k)
#print d[:,0]-k


P_window=np.array([.2,.2])  
C_window=.65	
nu=-2; n_pad=1000
# initialize the FASTPT class	

fastpt=FASTPT.FASTPT(k,nu,low_extrap=-5,high_extrap=5,n_pad=n_pad,verbose=True) 
		
t1=time()	
P_spt=fastpt.one_loop(P,C_window=C_window) 

t2=time()
print('time'), t2-t1 


print('To make a one-loop power spectrum for ', k.size, ' grid points, using FAST-PT takes ', t2-t1, 'seconds.')

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

fig=plt.figure(figsize=(16,10))

x1=10**(-2.5)
x2=10
ax1=fig.add_subplot(211)
ax1.set_ylim(1e-2,1e3)
ax1.set_xlim(x1,x2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$P_{22}(k)+ P_{13}(k)$ [Mpc/$h$]$^3$', size=30)
ax1.tick_params(axis='both', which='major', labelsize=30)
ax1.tick_params(axis='both', width=2, length=10)
ax1.tick_params(axis='both', which='minor', width=1, length=5)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax1.xaxis.labelpad = 20
ax1.set_xticklabels([])

ax1.plot(k,P_spt, lw=2,color='black', label=r'$P_{22}(k) + P_{13}(k)$, FAST-PT ' )
ax1.plot(k,-P_spt, '--',lw=2, color='black', alpha=.5 )


plt.grid()

ax2=fig.add_subplot(212)
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ [$h$/Mpc]', size=30)
ax2.set_ylim(.99,1.01)
ax2.set_xlim(x1,x2)
ax2.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', width=2, length=10)
ax2.tick_params(axis='both', which='minor', width=1, length=5)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax2.xaxis.labelpad = 20


ax2.plot(d[:,0],P_spt/(d[:,2]+d[:,3]),lw=2, color='black', alpha=.5, label='ratio to conventional method')

plt.legend(loc=3,fontsize=30)
plt.grid()

plt.tight_layout()
plt.show()
fig.savefig('example_plot.pdf')