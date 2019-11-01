''' This is the code I use to make animaitons of the RG flow. 
	The code reads in data in the same format RG_<integration type>.py 
	outputs. 
	
	You may get some warnings for taking log of zero and such. So, far these
	have not ruined any of the output 
	J. E. McEwen (c) 2016
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys 


# load data here
#name='data file to load'
#save_name=' animation to save , must be a .mp4 extension'
#title=r'what ever you like '
name='RG_STS_001_kmax10_ds1.npy'
save_name='test.mp4'
title=r'RG test'
d=np.load(name)


k=d[0,:]
P_0=d[2,:]


# downsample if number of frames is large
downsample=1
d=d[2::downsample,:]
N=d.shape[0]

# start program 
#Lambda=0

def n_eff(k,P):
	ln_P=np.log(P)
	ln_k=np.log(k) 
	return k[:-1], np.diff(ln_P)/np.diff(ln_k)

fig=plt.figure(figsize=(14,10))
fig.suptitle(title, size=20)
ax1=plt.subplot(1,2,1)
ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_ylim(80,1e5)
ax1.set_xlabel(r'$k$', size=20)
ax1.set_ylabel(r'$P_\star(k,\lambda)$', size=20)

ax1.plot(k,P_0, color='black', label=r'$P_\star(\lambda=0)$')
lambda_text=ax1.text(.05, .5, '', transform=ax1.transAxes, fontsize=40)

line1, =ax1.plot([], [], color='red', label=r'$P_\star(\lambda + \delta \lambda)$') 

plt.grid()			
plt.legend(loc=3)

ax2=plt.subplot(1,2,2)
ax2.set_xscale('log')
ax2.set_ylim(-3,.5)
ax2.set_xlabel(r'$k$', size=20)
ax2.set_ylabel(r'$n_{eff}(k)$', size=20)

a,b=n_eff(k,P_0)
ax2.plot(a,b, color='black', label=r'$P_\star(\lambda=0)$')


line2, =ax2.plot([], [], color='red') 
plt.axhline(-1.4, color='black', alpha=.5)

plt.grid()			


def init():
	line1.set_data([],[]) 
	line2.set_data([],[]) 
		
	return line1,line2

def animate(i):
	
	line1.set_data(k,d[i,:])
	a,b=n_eff(k,d[i,:])
	line2.set_data(a,b)
	
	Lambda=d[i,0]
	lambda_text.set_text(r'$\lambda$ = %2.4f' % Lambda)

	return line1,line2
	
ani=animation.FuncAnimation(fig, animate, init_func=init, frames=N, interval=1, blit=False) 
# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

# comment this out if you don't have ffmpeg 
ani.save(save_name, writer="ffmpeg")
 
	
			
plt.show() 