from cosmosis.datablock import names, option_section
import cPickle as pickle
import numpy as np
import string, pdb
from scipy import stats

# Many more redshift bin combinations than bins
# Get the lbin limits for each bin
# then apply the cuts.

def get_angular_frequency_cuts(block, samples, method_sh, method_pos, cut_per_bin, corr, chi_of_z, filename):
	lmin_ee, lmax_ee = [],[]
	lmin_nn, lmax_nn = [],[]
	lmin_ne, lmax_ne = [],[]

	if filename!=None:
		cuts_file = [np.loadtxt(filename).T[1], np.loadtxt(filename).T[2]]
	else:
		cuts_file = None

	if corr['ee']:
		nzbin = block[samples['shear'], 'nzbin']
		for i in xrange(nzbin):
			if cut_per_bin:
				lmin, lmax = choose_l_limits(block, cuts_file, i+1, samples['shear'], method_sh, chi_of_z )
			else:
				lmin, lmax = choose_l_limits(block, cuts_file, 1, samples['shear'], method_sh, chi_of_z)
			lmin_ee.append(lmin)
			lmax_ee.append(lmax)			

	if corr['nn']:
		nzbin = block[samples['pos'], 'nzbin']
		for i in xrange(nzbin):
			if cut_per_bin:
				lmin, lmax = choose_l_limits(block, cuts_file, i+1, samples['pos'], method_pos, chi_of_z )
			else:
				lmin,lmax = choose_l_limits(block, cuts_file, 1, samples['pos'], method_pos, chi_of_z )
			lmin_nn.append(lmin)
			lmax_nn.append(lmax)
	if corr['ne']:
		nzbin_1 = block[samples['pos'], 'nzbin']
		nzbin_2 = block[samples['shear'], 'nzbin']
		nzbin = max(nzbin_1, nzbin_2)
		for i in xrange(nzbin):
			if cut_per_bin:
				lmin, lmax = choose_l_limits(block, cuts_file, i+1, samples['pos'], method_pos, chi_of_z )
			else:
				lmin,lmax = choose_l_limits(block, cuts_file, 1, samples['pos'], method_pos, chi_of_z )
			lmin_ne.append(lmin)
			lmax_ne.append(lmax)

	return lmin_ee, lmax_ee,lmin_nn, lmax_nn, lmin_ne, lmax_ne

def apply_angular_frequency_cuts(block, samples, corr, lmin_ee, lmax_ee, lmin_nn, lmax_nn, lmin_ne, lmax_ne):

	correlations=[]
	if corr['ee']:
		nzbin_shear = block[samples['shear'], 'nzbin']
		correlations += [(nzbin_shear,nzbin_shear, 'galaxy_shape_cl', lmin_ee, lmax_ee)]	
	if corr['nn']:
		nzbin_pos = block[samples['pos'], 'nzbin']
		correlations += [(nzbin_pos,nzbin_pos, 'galaxy_position_cl', lmin_nn, lmax_nn)]
	if corr['ne']:
		nzbin_pos = block[samples['pos'], 'nzbin']
		nzbin_shear = block[samples['shear'], 'nzbin']
		correlations += [(nzbin_pos,nzbin_shear, 'galaxy_position_shape_cross_cl', lmin_ne, lmax_ne)]

	for  sp in correlations:
			maxbin1, maxbin2 = sp[0],sp[1]
			l = block[sp[2], 'ell']
			lmin, lmax = sp[3], sp[4]
			block[sp[2], 'ell_min']= lmin
			block[sp[2], 'ell_max']= lmax
			for i in xrange(maxbin1):	
				for j in xrange(maxbin2):
					
					upper =  lmax[j]
					lower =  lmin[j]
					print 'Applying scale cuts %e > ell > %e on bin %d %d'%(lower,upper,i+1,j+1)
					block[sp[2], 'lmin_%d_%d'%(i+1,j+1)]= lower
					block[sp[2], 'lmax_%d_%d'%(i+1,j+1)]= upper
					
def get_median_redshift(block, bin, cat):
	n_of_z = block[cat, 'bin_%d'%bin] 
	z = block[cat, 'z'] 
	prob_dist = n_of_z/n_of_z.sum()
	gen = stats.rv_discrete(values=(z,prob_dist), inc=z[1]-z[0])
	
	return gen.median()

def choose_l_limits(block, cuts_file, bin, sample, method, chi_of_z):
	if method=='fixed':
		return block[sample, 'lmin_%d'%bin], block[sample, 'lmax_%d'%bin]
	elif method=='rassat08':
		h = block['cosmological_parameters', 'h0']
		z_med = get_median_redshift(block, bin, sample)
		x = chi_of_z(z_med)

		kmax = 0.132 * z_med * h
		lmax = kmax * x
		return 10., lmax
	elif method=='file':
		return cuts_file[0][bin-1], cuts_file[1][bin-1]
