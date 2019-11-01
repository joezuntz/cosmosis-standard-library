from __future__ import division 
import numpy as np 
from matter_power_spt import one_loop
import FASTPT 
from time import time 

# load the input power spectrum data 

d_extend=np.loadtxt('inputs/P_lin.dat')
d=np.loadtxt('inputs/P_RSD_A.dat')
d2=np.loadtxt('inputs/P_RSD_B.dat')
k=d2[:,0]
P=d2[:,1]

k=d_extend[:-1,0]
P=d_extend[:-1,1]


# use if you want to interpolate data 
#from scipy.interpolate import interp1d 
#power=interp1d(k,P)
#k=np.logspace(np.log10(k[0]),np.log10(k[-1]),3000)
#P=power(k)
#print d[:,0]-k


P_window=np.array([.2,.2])  
C_window=.65	
n_pad=1000
# initialize the FASTPT class		
fastpt=FASTPT.FASTPT(k,to_do=['RSD'],low_extrap=-6,high_extrap=4,n_pad=n_pad) 
	
	
t1=time()	
f=1.
mu_n_list=[0.05,0.5,0.9]
mu_n=mu_n_list[0]
ABsum=fastpt.RSD_ABsum_mu(P,f,mu_n=mu_n,C_window=C_window) 
t2=time()

ABsum=ABsum[98:599]
print('To make a one-loop power spectrum for RSD', k.size, ' grid points, using FAST-PT takes ', t2-t1, 'seconds.')


from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times','Palatino']})
rc('text', usetex=True)


fig=plt.figure(figsize=(16,10))

x1=10**(-2.5)
x2=10
ax1=fig.add_subplot(231)
ax1.set_ylim(1e-2,1e4)
ax1.set_xlim(x1,x2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$A+B$ [Mpc/$h$]$^3$', size=25)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.tick_params(axis='both', width=2, length=10)
ax1.tick_params(axis='both', which='minor', width=1, length=5)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax1.xaxis.labelpad = 20
ax1.set_xticklabels([])

ax1.plot(k[98:599],ABsum, lw=3, color='black', label=r'$A+B$')
ax1.plot(k[98:599],-ABsum,'--', lw=3, color='black', label=r'$-(A+B)$')
ax1.set_title(r'$f=%2.1f$, $\mu_n=%2.2f$'%(f,mu_n),fontsize=20)

plt.legend(loc=2, fontsize=25)
plt.grid()

ax2=fig.add_subplot(234)
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ [$h$/Mpc]', size=25)
ax2.set_ylim(-.001,0.001)
ax2.set_xlim(x1,x2)

# labels = [item.get_text() for item in ax2.get_yticklabels()]
# labels[0] = r'$-1\times 10^{-3}$'
# labels[1] = r'$-5\times 10^{-4}$'
# labels[2] = '0'
# labels[3] = r'$5\times 10^{-4}$'
# labels[4] = r'$1\times 10^{-3}$'
# ax2.set_yticklabels(labels)

ax2.tick_params(axis='both', which='major', labelsize=25)
ax2.tick_params(axis='both', width=2, length=10)
ax2.tick_params(axis='both', which='minor', width=1, length=5)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax2.xaxis.labelpad = 20

A1=d[:,7]+f*d[:,8]
A3=f*d[:,9]+f**2 *d[:,10]
A5=f**2 *d[:,11]
Ap1=(d[:,2]+f*d[:,3])*d[:,1]
Ap3=(f*d[:,4]+f**2 *d[:,5])*d[:,1]
Ap5=f**2 *d[:,6]*d[:,1]
B0=d2[:,2]+2*f*d2[:,4]+f**2*d2[:,7]
B2=d2[:,3]+2*f*d2[:,5]+f**2*d2[:,8]
B4=2*f*d2[:,6]+f**2*d2[:,9]
B6=f**2*d2[:,10]
ABsum_mu2 = k[98:599]*f*(A1+Ap1) + (f*k[98:599])**2 *B0
ABsum_mu4 = k[98:599]*f*(A3+Ap3) + (f*k[98:599])**2 *B2
ABsum_mu6 = k[98:599]*f*(A5+Ap5) + (f*k[98:599])**2 *B4
ABsum_mu8 = (f*k[98:599])**2 *B6
jb_ABsum = ABsum_mu2 * mu_n**2 + ABsum_mu4*mu_n**4 + ABsum_mu6*mu_n**6 + ABsum_mu8*mu_n**8

ax2.plot(k[98:599],ABsum/jb_ABsum -1 ,lw=2, color='black')
ax2.text(0.07, 0.85, 'fractional difference',transform=ax2.transAxes,verticalalignment='bottom', horizontalalignment='left', fontsize=25, bbox=dict(facecolor='white', edgecolor='black', pad=8.0))

# plt.legend(loc=3,fontsize=20)
plt.grid()


mu_n=mu_n_list[1]
ABsum=fastpt.RSD_ABsum_mu(P,f,mu_n=mu_n,C_window=C_window) 
ABsum=ABsum[98:599]
ax3=fig.add_subplot(232)
ax3.set_ylim(1e-2,1e4)
ax3.set_xlim(x1,x2)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.tick_params(axis='both', which='major', labelsize=25)
ax3.tick_params(axis='both', width=2, length=10)
ax3.tick_params(axis='both', which='minor', width=1, length=5)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax3.xaxis.labelpad = 20
ax3.set_xticklabels([])
ax3.set_yticklabels([])

ax3.plot(k[98:599],ABsum, lw=3, color='black')
ax3.plot(k[98:599],-ABsum,'--', lw=3, color='black')
ax3.set_title(r'$f=%2.1f$, $\mu_n=%2.2f$'%(f,mu_n),fontsize=20)
# plt.legend(loc=3, fontsize=25)
plt.grid()

ax4=fig.add_subplot(235)
ax4.set_xscale('log')
ax4.set_xlabel(r'$k$ [$h$/Mpc]', size=25)
# ax4.set_ylabel('fractional difference',size=30)
ax4.set_ylim(-.001,0.001)
ax4.set_xlim(x1,x2)
ax4.tick_params(axis='both', which='major', labelsize=25)
ax4.tick_params(axis='both', width=2, length=10)
ax4.tick_params(axis='both', which='minor', width=1, length=5)
ax4.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax4.xaxis.labelpad = 20
ax4.set_yticklabels([])

jb_ABsum = ABsum_mu2 * mu_n**2 + ABsum_mu4*mu_n**4 + ABsum_mu6*mu_n**6 + ABsum_mu8*mu_n**8

ax4.plot(k[98:599],ABsum/jb_ABsum -1 ,lw=2, color='black')

# plt.legend(loc=3,fontsize=20)
plt.grid()

mu_n=mu_n_list[2]
ABsum=fastpt.RSD_ABsum_mu(P,f,mu_n=mu_n,C_window=C_window) 
ABsum=ABsum[98:599]

ax5=fig.add_subplot(233)
ax5.set_ylim(1e-2,1e4)
ax5.set_xlim(x1,x2)
ax5.set_xscale('log')
ax5.set_yscale('log')
# ax5.set_ylabel(r'RSD $A+B$', size=30)
ax5.tick_params(axis='both', which='major', labelsize=25)
ax5.tick_params(axis='both', width=2, length=10)
ax5.tick_params(axis='both', which='minor', width=1, length=5)
ax5.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax5.xaxis.labelpad = 20
ax5.set_xticklabels([])
ax5.set_yticklabels([])

ax5.plot(k[98:599],ABsum, lw=3, color='black')
ax5.plot(k[98:599],-ABsum,'--', lw=3, color='black')
ax5.set_title(r'$f=%2.1f$, $\mu_n=%2.2f$'%(f,mu_n),fontsize=20)
# plt.legend(loc=3, fontsize=25)
plt.grid()

ax6=fig.add_subplot(236)
ax6.set_xscale('log')
ax6.set_xlabel(r'$k$ [$h$/Mpc]', size=25)
# ax6.set_ylabel('fractional difference',size=30)
ax6.set_ylim(-.001,0.001)
ax6.set_xlim(x1,x2)
ax6.tick_params(axis='both', which='major', labelsize=25)
ax6.tick_params(axis='both', width=2, length=10)
ax6.tick_params(axis='both', which='minor', width=1, length=5)
ax6.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax6.xaxis.labelpad = 20
ax6.set_yticklabels([])

jb_ABsum = ABsum_mu2 * mu_n**2 + ABsum_mu4*mu_n**4 + ABsum_mu6*mu_n**6 + ABsum_mu8*mu_n**8

ax6.plot(k[98:599],ABsum/jb_ABsum -1 ,lw=2, color='black')

# plt.legend(loc=3,fontsize=20)
plt.grid()

plt.tight_layout()
plt.show()
fig.savefig('RSD_ABsum_plot.pdf')