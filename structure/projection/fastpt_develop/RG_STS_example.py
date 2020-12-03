''' 
    Example script to run the STS method for RG results. 
    We compare these results to the same resulst from the Cotper code.
'''

import numpy as np
from RG_STS import RG_STS

# load the Copter data 
d=np.loadtxt('RGcopter_10_1000.dat')
k=d[:,0]; P=d[:,1]; copter=d[:,2]

# set the STS parameters here 
# this combo seems to work well for k_max=10, 2000 grid points 
# if you encounter instabilities you may want to fiddle with these values. 
stage=10
mu=.1
# this is \Delta \lambda_{CFL}
dlambda_CFL=1e-3

STS_params=[stage, mu, dlambda_CFL] 

P_window=np.array([.2,.2]) 
C_window=.75 
step=.01
max=1
n_pad=500
P_rg=RG_STS(STS_params,'test_STS',k,P,step,max,n_pad,P_window,C_window)	

import matplotlib.pyplot as plt

ax=plt.subplot(121)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$k$', size=20)
ax.set_ylabel(r'$P(k)$', size=20) 

ax.plot(k,P, label='linear')
ax.plot(k,copter, label='RG copter')
ax.plot(k, P_rg, label='RG FAST-PT')

plt.grid()
plt.legend(loc=3)

ax=plt.subplot(122)
ax.set_xscale('log')
ax.set_ylim(.8,1.2)
ax.set_xlabel(r'$k$', size=20)

ax.plot(k,P_rg/copter, label='ratio')

plt.grid()
plt.legend(loc=3)

plt.show()