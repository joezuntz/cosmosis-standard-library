''' 
    Example script to plot the 
    galaxy componets using the FASTPT algorithm.
    
    Additionally the code also makes use of extending the power spectrum out
    to higher and lower k to alliviate edge effects. 
    
    The parameters low_extrap and high_extrap control the extension of P_{lin}(k) 
    out to lower and higher k values by extrapolating using the effective power-law
    index n_{eff}= dlogP/dlogk. 
''' 


import numpy as np
# example data 
d=np.loadtxt('P_bias_example.dat')
#print('shape', d.shape)
'''
0: k [in h/Mpc]
1: Plin (scales as D^2)
2: P22 (scales as D^4)
3: P13 (scales as D^4)
4: b1 b2 (scales as D^4)
5: b2^2 (scales as D^4)
6: bs^2 (scales as D^4)
7: b1 bv (no advection) (scales as D^2)
8: b1 bv Ls (scales as D^2)
9: b2 bv (scales as D^2)
10: bs bv (scales as D^2)
11: bv^2 (constant )
'''

# example use of FASTPTII
k=d[:,0]; P_lin=d[:,1]

import FASTPT
from time import time 
C_window=.65	
n_pad=1000
		
t1=time()	
# initialize the FASTPT class
nu=-2
fastpt=FASTPT.FASTPT(k,-2,n_pad=n_pad) 

P_1loop=fastpt.one_loop(P_lin,C_window=C_window) 
_,Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig4 =fastpt.P_bias(P_lin,C_window=C_window) 
t2=time()
print('The time to make density-density type power spectra is ', t2-t1,' .')


# plot the output 
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
gs = gridspec.GridSpec(2,3, height_ratios=[2,1.25])


fig=plt.figure(figsize=(16,10))

#x1=10**(-2.5)
#x2=50
x1=10**(-5)
x2=100
x1=k[0]
x2=k[-1]
ax1=fig.add_subplot(gs[0,0])
#ax1.set_ylim(1e-2,1e3)
ax1.set_xlim(x1,x2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$P(k)$ [Mpc/$h$]$^3$', size=20)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', width=2, length=10)
ax1.tick_params(axis='both', which='minor', width=1, length=5)
#ax1.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax1.xaxis.labelpad = 20
ax1.set_xticklabels([])

ax1.plot(k,P_1loop, lw=2,color='black', label=r'$P_{22}(k) + P_{13}(k)$, FAST-PT ' )
ax1.plot(k,-P_1loop, '--',lw=2, color='black', alpha=.5 )
plt.legend(loc=3)
plt.grid()

ax2=fig.add_subplot(gs[1,0])
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ [$h$/Mpc]', size=20)
ax2.set_ylabel('ratio', size=20)
ax2.set_ylim(.99,1.01)
ax2.set_xlim(x1,x2)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', width=2, length=10)
ax2.tick_params(axis='both', which='minor', width=1, length=5)
#ax2.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax2.xaxis.labelpad = 20
 
 
ax2.plot(d[:,0],P_1loop/(d[:,2]+d[:,3]),lw=2, color='black', alpha=.5)
plt.grid()

##########################################################################
ax1=fig.add_subplot(gs[0,1])
#ax1.set_ylim(1e-2,1e3)
ax1.set_xlim(x1,x2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$P(k)$ [Mpc/$h$]$^3$', size=20)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', width=2, length=10)
ax1.tick_params(axis='both', which='minor', width=1, length=5)
#ax1.xaxis.set_major_formatter(FormatStrFormatter('%2.4f'))
ax1.xaxis.labelpad = 20
ax1.set_xticklabels([])

ax1.plot(k,Pd1d2, lw=2, color='black',label=r'$b_1 b_2$, FAST-PT ' )
ax1.plot(k,Pd2d2, '--', lw=2,color='red', label=r'$b_2^2$, FAST-PT ' )
plt.legend(loc=3)

plt.grid()

ax2=fig.add_subplot(gs[1,1])
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ [$h$/Mpc]', size=20)
ax2.set_ylabel('ratio', size=20)
ax2.set_ylim(.99,1.01)
ax2.set_xlim(x1,x2)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', width=2, length=10)
ax2.tick_params(axis='both', which='minor', width=1, length=5)
#ax2.xaxis.set_major_formatter(FormatStrFormatter('%2.4f'))
ax2.xaxis.labelpad = 20
 
 
ax2.plot(k,Pd1d2/d[:,4], lw=2,color='black', alpha=.5, label=r'$b_1 b_2$, FAST-PT ' )
ax2.plot(k,Pd2d2/d[:,5]/4., '--', lw=2,color='red', label=r'$b_2^2$, FAST-PT ' )
plt.grid()

##########################################################################
ax1=fig.add_subplot(gs[0,2])
#ax1.set_ylim(1e-2,1e3)
ax1.set_xlim(x1,x2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$P(k)$ [Mpc/$h$]$^3$', size=20)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', width=2, length=10)
ax1.tick_params(axis='both', which='minor', width=1, length=5)
#ax1.xaxis.set_major_formatter(FormatStrFormatter('%2.4f'))
ax1.xaxis.labelpad = 20
ax1.set_xticklabels([])

ax1.plot(k,-Pd1s2, lw=2, color='black', label=r'$-b_1 b_s$, FAST-PT ' )
ax1.plot(k,Pd2s2, '--', lw=2,color='red', label=r'$b_2 b_s$, FAST-PT ' )
ax1.plot(k,Ps2s2, '-.', lw=2,color='blue', label=r'$b_s^2$, FAST-PT ' )
plt.legend(loc=3)

plt.grid()

ax2=fig.add_subplot(gs[1,2])
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ [$h$/Mpc]', size=20)
ax2.set_ylabel('ratio', size=20)
ax2.set_ylim(.99,1.01)
ax2.set_xlim(x1,x2)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', width=2, length=10)
ax2.tick_params(axis='both', which='minor', width=1, length=5)
#ax2.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax2.xaxis.labelpad = 20
 
ax2.plot(k,Pd1s2/d[:,12], lw=2, color='black', alpha=.5, label=r'$-b_1 b_s$, FAST-PT ' )
ax2.plot(k,Pd2s2/d[:,13]/2., '--', lw=2,color='red', label=r'$b_2 b_s$, FAST-PT ' )
ax2.plot(k,Ps2s2/d[:,6]/4., '-.', lw=2,color='blue', label=r'$b_s^2$, FAST-PT ' )
plt.grid()

plt.tight_layout()
plt.show()
fig.savefig('bias_example_plot.pdf')