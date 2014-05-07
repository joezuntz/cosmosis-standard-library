import os
import pydesglue
import numpy as np
try:
	import pylab
except ImportError:
	pylab=None


class Rosenbrock(object):
	def __init__(self):
	        self.a1 = 100.0
        	self.a2 = 20.0

	def logl(self,p):
		return -(self.a1 * (p[1] - p[0] ** 2) ** 2 + (1 - p[0]) ** 2) / self.a2

	def __call__(self, p):
		return logl(self,p)
		


def execute(handle):
	try:
		package = pydesglue.DesDataPackage.from_fits_handle(handle)
		section = pydesglue.section_names.cosmological_parameters
		p=[]
		p.append(package.get_param(section,"P1"))
		p.append(package.get_param(section,"P2"))

		rosenbrock = Rosenbrock()
		#compute likelihood
		like =rosenbrock.logl(p)
	
		section = pydesglue.section_names.likelihoods
		package.set_param(section,"ROSENBROCK_LIKE",like)
		package.write_to_fits_handle(handle)
	
	except KeyboardInterrupt:
		#If the user presses Ctrl-C we respect that
		raise KeyboardInterrupt
	except Exception as E:
		#But any other kind of error is a problem
                print "There was an error calculating the Rosenbrock likelihood: "
                return 1
        return 0

