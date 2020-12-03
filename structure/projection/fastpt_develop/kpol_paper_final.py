from __future__ import division
import numpy as np 
from matter_power_spt import one_loop
import FASTPT
from time import time 

# load the input power spectrum data 
d=np.genfromtxt('inputs/P_OV_KPol.dat',skip_header=1)
k=d[:,0]
P=d[:,1]


# d_extend=np.loadtxt('PT_Plininterp4_z0.dat')
d_extend=np.genfromtxt('inputs/P_lin.dat',skip_header=1)
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
fastpt=FASTPT.FASTPT(k,to_do=['kPol'],low_extrap=-6,high_extrap=4,n_pad=n_pad) 
	
	
t1=time()	
p1,p2,p3=fastpt.kPol(P,C_window=C_window) 
t2=time()

k=k[98:599]
p1,p2,p3=p1[98:599],p2[98:599],p3[98:599]


print('To make a one-loop power spectrum for ', k.size, ' grid points, using FAST-PT takes ', t2-t1, 'seconds.')

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
ax1.set_ylim(1e-8,1e8)
ax1.set_xlim(x1,x2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('$P^{(m)}(k)$ [Mpc/$h$]$^7$', size=25)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.tick_params(axis='both', width=2, length=10)
ax1.tick_params(axis='both', which='minor', width=1, length=5)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax1.xaxis.labelpad = 20
ax1.set_xticklabels([])

ax1.plot(k,p1, lw=4, color='black', label=r'$P^{(0)}(k)$')
plt.legend(loc=3,fontsize=25)
plt.grid()


y1=-0.00007
y2=-y1
ax2=fig.add_subplot(234)
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ [$h$/Mpc]', size=25)
ax2.set_ylim(y1,y2)
ax2.set_xlim(x1,x2)

# labels = [item.get_text() for item in ax2.get_yticklabels()]
# # labels[0] = r'$-8\times 10^{-5}$'
# labels[1] = r'$-6\times 10^{-5}$'
# labels[2] = r'$-4\times 10^{-5}$'
# labels[3] = r'$-2\times 10^{-5}$'
# labels[4] = '0'
# labels[5] = r'$2\times 10^{-5}$'
# labels[6] = r'$4\times 10^{-5}$'
# labels[7] = r'$6\times 10^{-5}$'
# # labels[8] = r'$8\times 10^{-5}$'
# ax2.set_yticklabels(labels)

ax2.tick_params(axis='both', which='major', labelsize=25)
ax2.tick_params(axis='both', width=2, length=10)
ax2.tick_params(axis='both', which='minor', width=1, length=5)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax2.xaxis.labelpad = 20

ax2.plot(k,p1/d[:,5]-1,lw=2, color='black')
ax2.text(0.07, 0.07, 'fractional difference',transform=ax2.transAxes,verticalalignment='bottom', horizontalalignment='left', fontsize=25, bbox=dict(facecolor='white', edgecolor='black', pad=8.0))
# plt.legend(loc=3,fontsize=20)
plt.grid()

ax3=fig.add_subplot(232)
ax3.set_ylim(1e-8,1e8)
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

ax3.plot(k,p2, lw=4, color='black', label=r'$P^{(\pm 1)}(k)$')
# plt.legend(loc=3,fontsize=25)
plt.grid()

ax4=fig.add_subplot(235)
ax4.set_xscale('log')
ax4.set_xlabel(r'$k$ [$h$/Mpc]', size=25)
ax4.set_ylim(y1,y2)
ax4.set_xlim(x1,x2)
ax4.tick_params(axis='both', which='major', labelsize=25)
ax4.tick_params(axis='both', width=2, length=10)
ax4.tick_params(axis='both', which='minor', width=1, length=5)
ax4.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax4.xaxis.labelpad = 20
ax4.set_yticklabels([])

ax4.plot(k,p2/d[:,7]-1,lw=2, color='black')

# plt.legend(loc=3,fontsize=20)
plt.grid()

ax5=fig.add_subplot(233)
ax5.set_ylim(1e-8,1e8)
ax5.set_xlim(x1,x2)
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.tick_params(axis='both', which='major', labelsize=25)
ax5.tick_params(axis='both', width=2, length=10)
ax5.tick_params(axis='both', which='minor', width=1, length=5)
ax5.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax5.xaxis.labelpad = 20
ax5.set_xticklabels([])
ax5.set_yticklabels([])

ax5.plot(k,p3, lw=4, color='black', label=r'$P^{(\pm 2)}(k)$')
# plt.legend(loc=3,fontsize=25)
plt.grid()

ax6=fig.add_subplot(236)
ax6.set_xscale('log')
ax6.set_xlabel(r'$k$ [$h$/Mpc]', size=25)
ax6.set_ylim(y1,y2)
ax6.set_xlim(x1,x2)
ax6.tick_params(axis='both', which='major', labelsize=25)
ax6.tick_params(axis='both', width=2, length=10)
ax6.tick_params(axis='both', which='minor', width=1, length=5)
ax6.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
ax6.xaxis.labelpad = 20
ax6.set_yticklabels([])

ax6.plot(k,p3/d[:,9]-1,lw=2, color='black')

# plt.legend(loc=3,fontsize=20)
plt.grid()

plt.tight_layout()
plt.show()
fig.savefig('kpol_plot.pdf')