#Code by Donnacha Kirk
#Edited by Simon Samuroff 09/2015 

import numpy as np
import pdb, matplotlib.pyplot as plt
import cPickle as pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import pdb

def setup(options):
	dz = options.get_double(option_section, "dz", default=0.01)
	survey = options.get_string(option_section, "survey")
	try: 
		cat_opt = options.get_string(option_section, "catastrophic_outliers")
		cat_opt = {'mode': cat_opt, 'fcat': block[survey,'fcat'], 'dzcat': block[survey,'dzcat'],  'zcat0': block[survey,'zcat0'],  'zcat': block[survey,'zcat'],  'sigcat': block[survey,'sigcat'] }
	except: cat_opt = None
	
	return (dz, survey, cat_opt)

def gaussian(z,mu,sigma):
	return np.exp(-0.5*(z-mu)**2/sigma**2) / np.sqrt(2*np.pi) / sigma

def smail_distribution(z, alpha, beta, z0):
	return (z**alpha) * np.exp(-(z/z0)**beta)

def sigma_phot(sigma_z, alpha_z, z):
	return sigma_z*(1+z)**alpha_z

def photometric_error(z, Nz, sigma_z, alpha_z, bias, cat_opt):
	# Calculates the n(zphot|zspec), assuming a Gaussian error distribution
	# Also includes catastrophic outliers
	# See Hearin et al (2010) for the prescription used
	nz = len(z)
	output = np.zeros((nz,nz))

	cat_mode = None
	if cat_opt!=None:
		cat_mode = cat_opt['mode']
		fcat = cat_opt['fcat']
		dzcat = cat_opt['dzcat']
		zcat = cat_opt['zcat']
		sigcat= cat_opt['sigcat']
		try: zcat0= cat_opt['zcat0']
		except:
			zcat0 = np.random.random()*(z[-1]-z[0])+z[0]

	dz = (z[1]-z[0])

	Nz /= sum(Nz * (z[1]-z[0]))

	for i in xrange(nz):
		p = gaussian(z,z[i]-bias, sigma_phot(sigma_z, alpha_z, z[i]) )
		
		# Subtract some probability from a given region of the p(z)
		if cat_opt!=None:
			step=(dzcat/2.)-abs(z-zcat0)
			step[step==0.0]=-1.0
			step=0.5*(step/abs(step)+1.0)

		# Define a Gaussian island of outliers, normalised to 
		# the probability scattered from the affected region
		if cat_mode=='island':
			pcat = ((2.0*np.pi)**-0.5 / sigcat) * np.exp(-1.0*(z-zcat)**2 / (2.*sigcat*sigcat))
			pcat = pcat * sum(step*fcat*p*dz) / sum(pcat*dz)
			p = (1.-step*fcat)*p + pcat

		# Or scatter it uniformly across all z
		if cat_mode=='uniform':
			p=(1.-step*fcat)*p + step*fcat/(z[-1]-z[0])
		
		# Normalise
		p /= sum(p*(z[1]-z[0]))

		output[:,i] = (p * Nz[i])

	w = np.where((z>0.) & (z<0.511))[0]
	return output

def find_bins(z, nz_true, nbin):
	nz_true = nz_true/nz_true.sum()*nbin
	cum = np.cumsum(nz_true)
	bin_edges = [0.0]
	for i in xrange(1,nbin):
		edge = np.interp(1.0*i,cum,z)
		bin_edges.append(edge)
	bin_edges.append(z.max())	
	return np.array(bin_edges)

def compute_bin_nz(z_prob_matrix, z, edges, ngal):
	NI = []
	nbin = len(edges)-1
	dz = z[1]-z[0]
	for low,high in zip(edges[:-1], edges[1:]):
		w = np.where((z>low) & (z<high))[0]
		# Sum over all possible ztrue
		ni = z_prob_matrix[w,:].sum(axis=0)

		# Normalise the n(z) in each redshift bin to 1 over the redshift range
		# of the survey
		ni*= 1.0/(ni.sum()*dz)							#ngal/(nbin*ni.sum()*dz)
		assert(len(ni)==len(z))
		NI.append(ni)
	return NI
	
def compute_nz(smail_dist, z, nbin, sigma_z, alpha_z, ngal, bias, cat_opt):
	#Set up Smail distribution of z vector as the distribution of true redshifts of the galaxies, n(ztrue)
	nz_true = smail_dist #ribution(z, alpha, beta, z0)

	# Multiply that by a Gaussian to get the probability distribution of the measured photo-z for each true redshift
	# This gives a 2D probability distribution
	z_prob_matrix = photometric_error(z, nz_true, sigma_z, alpha_z, bias, cat_opt)
	edges = find_bins(z,nz_true,nbin)
	bin_nz = compute_bin_nz(z_prob_matrix, z, edges, ngal)
	return edges,bin_nz

def execute(block, config):
	dz, survey = config[0:2]
	# Keep the catastrophic outlier options together
	cat_opt = config[2]
	params = section_names.number_density_params
	alpha = block[survey, "smail_alpha"]
	beta = block[survey, "smail_beta"]
	z0 = block[survey, "z0"]

	sigma_z = block.get_double(survey, "sigz")
	alpha_z = block.get_double(survey, "alpha_z")
	bias = block.get(survey, "photoz_bias")
	ngal = block.get_double(survey, "ngal")
	zmax = block.get_double(survey, "zmax")
	nbin = block.get_int(survey, "nzbin")

	
	#Compute the redshift vector
	z = np.arange(0.0,zmax+dz/2,dz)
	
	dist = smail_distribution(z, alpha, beta, z0)
	#Run the main code for getting n(z) in bins
	edges,bins = compute_nz(dist, z, nbin, sigma_z, alpha_z, ngal, bias, cat_opt)

	#Save the results
	block[survey,"nz"] = len(z)
	block[survey,"z"] = z

	#Loop through the bins
	for i,bin in enumerate(bins):
		#The bin numbering starts at 1
		b=i+1
		name = "bin_%d" % b
		#Save the bin edges as parameters
		block[survey,"edge_%d"%b] = edges[i]
		#And save the bin n(z) as a column
		block[survey, name] =  bin
	#Also save the upper limit to the top bin
	block[survey,"edge_%d"%(nbin+1)] = edges[-1]

	return 0
		
def cleanup(config):
	#nothing to do here!  We just include this 
	# for completeness
	return 0
