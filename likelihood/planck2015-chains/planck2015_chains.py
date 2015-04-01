"""
At the time of writing the Planck collaboration has released 
chains for their 2015 data, but not the likelihood code that
generated them.

Until the latter is available we can use the released chains
to get an approximate and reasonable likelihood as long as
the assumptions and parameters made in the chains match what
the user wants to do (in terms of datasets, physics choices, 
etc).

This temporary module uses Kernel Density Estimation to 
approximate the likelihood from the chain.

It can be run with any of the released Planck chains
and any parameters they contains - you need to download some 
planck chains (they come in sets of 4) and the corresponding
.paramnames file, and set the "base" parameter to point
to the directory and base name they are in.

The default assumes you downloaded base_plikHM_TTTEEE_lowTEB 
files from the large fullGrid download to the directory this
file is in.

"""

from cosmosis.datablock import names, option_section
from scipy.stats.kde import gaussian_kde
import numpy as np
import os


cosmo=names.cosmological_parameters

def cosmosis_to_planck(block, names):
	mapping = {
		"omegabh2":"ombh2",
		"omegach2":"omch2",
		"ns":"n_s",
		"H0":"hubble",
		"sigma8":"sigma_8",
		"omegal":"omega_lambda",
		"omegam":"omega_m",
		"omegamh2":"ommh2",
		"omeganuh2":"omnuh2",
		"tau":"tau",
	}
	values = []
	for name in names:
		if name in mapping:
			values.append(block[cosmo,mapping[name]])
		elif name=='logA':
			values.append(np.log(10**10 * block[cosmo,'A_s']))
		else:
			raise ValueError("Do not know how to find Planck parameter called '%s' in cosmosis"%name)
	return values

class Temp(object):
	def get(self,section, name, default):
		return default

def setup(options):
	#look in current dir by default
	dirname = os.path.split(__file__)[0]+os.path.sep
	#Decide which file and columns in that file to use in KDE
	base = options.get_string(option_section, "base", default=dirname+"base_plikHM_TTTEEE_lowTEB")
	cols = options.get_string(option_section, "cols", default="omegabh2 omegach2 logA ns H0").split()
	#Look up names of columns in paramnames file that comes along with planck
	file_cols = [line.split()[0].strip('*') for line in open(base+'.paramnames')]
	#columns in chain are shifted by two places (for multiplicicty and likelihood)
	#from paramnames
	indices = [file_cols.index(c)+2 for c in cols]
	#read in chains, extracting only desired columns
	data_files = [os.path.join(dirname, "%s_%d.txt"%(base,i)) for i in [1,2,3,4]]
	data = [np.loadtxt(filename,usecols=indices).T for filename in data_files]
	#stack and build KDE
	data = np.hstack(data)
	kde = gaussian_kde(data)
	print
	print "The Planck 2015 chain KDE code has not been tested at all." 
	print
	return {'kde':kde, 'cols':cols}


def execute(block, config):
	kde = config['kde']
	cols = config['cols']
	x = cosmosis_to_planck(block, cols)
	like = np.log(kde.evaluate(x)[0])
	block[names.likelihoods, "planck2015_like"] = like
	return 0
