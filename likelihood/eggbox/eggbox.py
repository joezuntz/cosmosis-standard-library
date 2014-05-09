import os
import pydesglue
import numpy as np
try:
	import pylab
except ImportError:
	pylab=None


class Eggbox(object):
	"""Adapted from https://github.com/dfm/emcee/blob/master/examples/eggbox.py """
	def __init__(self):
		self.tmax = 10.0*np.pi
		self.constant = np.log(1.0/(self.tmax*self.tmax))

	def logprior(self,t):
		if (t[0] > self.tmax or t[0] < -self.tmax or \
		t[1] > self.tmax or t[1] < -self.tmax):
			return -np.inf
		else:
			return self.constant

	def loglhood(self,t):
		return (2.0 + np.cos(t[0]/2.0)*np.cos(t[1]/2.0))**5.0 + self.logprior(t)

	def __call__(self, t):
		return self.loglhood(t)

		


def execute(handle):
	try:
		package = pydesglue.DesDataPackage.from_fits_handle(handle)
		section = pydesglue.section_names.cosmological_parameters
		t=[]
		t.append(package.get_param(section,"T1"))
		t.append(package.get_param(section,"T2"))

		eggbox = Eggbox()
		#compute likelihood
		like =eggbox.loglhood(t)
	
		section = pydesglue.section_names.likelihoods
		package.set_param(section,"EGGBOX_LIKE",like)
		package.write_to_fits_handle(handle)
	
	except KeyboardInterrupt:
		#If the user presses Ctrl-C we respect that
		raise KeyboardInterrupt
	except Exception as E:
		#But any other kind of error is a problem
                print "There was an error calculating the Eggbox likelihood: "
                return 1
        return 0

