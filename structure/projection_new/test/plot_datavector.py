import pylab
import numpy as np
import sys
import os
import twopoint

spec_names={'wtheta': r'$w(\theta)$',
			'xip': r'$\xi_+$',
			'xim': r'$\xi_-$',
			'gammat': r'$\gamma_t$'}

tpt_orig=twopoint.TwoPointFile.from_fits('simulated_y1_fiducial_v1.fits',covmat_name=None)
tpt_new=twopoint.TwoPointFile.from_fits('simulated_y1_fiducial_v1_new.fits',covmat_name=None)

fig=pylab.figure()
ax=fig.add_subplot(111)

dv_ind_start=0
dv_starts=[]
for s in tpt_orig.spectra:
	l = len(s.value)
	ax.plot(np.arange(dv_ind_start,dv_ind_start+l), (tpt_new.get_spectrum(s.name).value)/s.value)
	dv_starts.append(dv_ind_start)
	dv_ind_start+=l

dv_length=dv_ind_start
ylims=ax.get_ylim()
for l in dv_starts:
	ax.plot([l,l], ylims, 'k-')
ax.set_ylim(ylims)
ax.set_xlim([0,dv_length])
ax.set_ylabel('new / orginal')
try:
	fig.savefig(sys.argv[1])
except IndexError:
	pylab.show()


    


