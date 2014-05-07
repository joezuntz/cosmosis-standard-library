import os
import sys
import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


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



HST_H0_MEAN = 0.738
HST_H0_SIGMA = 0.024

def setup(options):
	section = option_section
	mean = options.get_double(section, "mean", default=HST_H0_MEAN)
	sigma = options.get_double(section, "sigma", default=HST_H0_SIGMA)
	norm = 0.5
	return (mean, sigma, norm)


def execute(block, config):
	# Configuration data, read from ini file above
	t1 = block[cosmo, 'T1']
	t2 = block[cosmo, 'T2']
	print t1,t2
 	t=[t1,t2]
	eggbox = Eggbox()
	#compute likelihood
	like =eggbox.loglhood(t)
	block[likes, 'EGGBOX_LIKE'] = like


def cleanup(config):
	#nothing to do here!  We just include this
	# for completeness
	return 0
