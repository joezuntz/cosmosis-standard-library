from gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names, option_section
import cPickle as pickle
import numpy as np
import string
import pdb

# NOTE:
# The covaraiance matrix calculation requires some bookkeeping to keep track of
# 4 redshift bins i,j,k,m and two types of spectrum. 
# The full covariance matrix has the following structure:
# - nine blocks, each corresponding to the covariance for a pair of spectrum types
# 	- nzpair x nzpair sub-blocks, each corresponding to a matrix for a pair of bin combinations (i,j), (k,m)
#		- one nlbin x nlbin matrix
# build_covariance is called once and loops over spectrum type block
# compute_covariance is called once for each block and loops over redshift bin combinations
# Each matrix is written into place in a blank sub-block, which is then fitted into a block
# Finally the blocks are fitted into the full (nlbin x nzpair x nspect) x (nlbin x nzpair x nspect) covariance matrix

# Also note that the datavector does not contain equal numbers of the different types of spectra
# It contains nzbin*nzbin ne terms, but only nzbin*(nzbin+1)/2 of the ee and nn terms
# As such, looking for a specific block is not just a matter of splitting the matrix into three along each axis
# Rather, cov_eeee = cov[                : nl*nz*(nz+1)/2,                           : nl*nz*(nz+1)/2 ],
#	  cov_nnnn = cov[ nl*nz*(nz+1)/2 : nl*nz*(nz+1),              nl*nz*(nz+1)/2 : nl*nz*(nz+1)],
# 	  cov_enen = cov[ nl*nz*(nz+1)   : nl*nz*(nz+1) + nl*nz*nz,     nl*nz*(nz+1) : nl*nz*(nz+1) + nl*nz*nz] 



cl_sections= {'ee':'galaxy_shape_cl', 'ne':'galaxy_position_shape_cross_cl', 'nn':'galaxy_position_cl'}

class ClLikelihood(GaussianLikelihood):
	constant_covariance=False
	#one section name
	#don't define section name when reading options
	def __init__(self, options):

		constcov = options.get_bool(option_section, 'fixed_covariance')
		if constcov:
			self.constant_covariance=True
			covmat_file = options.get_string(option_section, 'covariance')
		data_file = options.get_string(option_section, 'data')
		self.shear_cat = options.get_string(option_section, 'shear_sample')
		self.pos_cat = options.get_string(option_section, 'clustering_sample')
		self.shear = options.get_bool(option_section, 'shear')
		self.pos = options.get_bool(option_section, 'galaxy_clustering')
		self.ggl = options.get_bool(option_section, 'ggl')
		auto = options.get_bool(option_section, 'auto_zbins')
		cross = options.get_bool(option_section, 'cross_zbins')

		self.interpolate = options.get_bool(option_section, 'interpolate')

		# Put the information specifying the datavector where
		# it can be accessed later
		self.spectra_to_use = []
		if self.shear:	self.spectra_to_use.append('ee')
		if self.ggl:	self.spectra_to_use.append('ne')
		if self.pos:	self.spectra_to_use.append('nn')

		self.x_section = [cl_sections[sp] for sp in self.spectra_to_use]
		self.x_name    = 'l_bin_edges'
		self.y_section = [cl_sections[sp] for sp in self.spectra_to_use]

		if len(self.spectra_to_use)==1:
			names.data_vector = "cl_%s"%self.spectra_to_use[0]
			self.like_name = "cl_%s"%self.spectra_to_use[0]
		else:
			names.data_vector = 'combined_cl_%s'%string.join(self.spectra_to_use).replace(' ','_')
			self.like_name = 'combined_cl_%s'%string.join(self.spectra_to_use).replace(' ','_')

		# Load the datavector and covariance matrix if fixed from the disc
		self.data_x, self.data_y = self.build_data(data_file)
		self.data_y_uncut = np.copy(self.data_y)
		if self.constant_covariance:
			self.covmat_file = covmat_file
			self.extract_covariance()

	def initialise_theory(self, block, auto=None, cross=None):
		"""
		Sets up the theory datavector and covariance matrix, if not using a
		fixed one. These are cosmology dependent and so could not be loaded
		in the setup function at the start of the pipeline.
		"""

		# Get the reshift binning
		self.get_zbin_combinations(block, auto, cross)

		# Choose which Cls to load
		self.spectra_theory = []
		self.n_elem = 0
		for sp in self.spectra_to_use:
			bins = getattr(self, 'bins_%s'%sp)
			for zcomb in bins: 
				# For the symmetric spectra ee and nn, only half the redshift
				# bin combinations are required
				if (sp[0]==sp[1]) and (zcomb[1]<zcomb[0]): continue
				else:
					self.spectra_theory.append((sp, zcomb))
					l = block[cl_sections[sp],'ell']
					if self.cuts:
						scale_cut = (l > block[cl_sections[sp],'lmin_%s_%s'%(zcomb[0]+1,zcomb[1]+1) ]) & (l < block[cl_sections[sp],'lmax_%s_%s'%(zcomb[0]+1,zcomb[1]+1)])
					else:
						scale_cut = np.ones_like(l).astype(np.bool)
					self.n_elem += len(l[scale_cut])
		# This list should now fully specify the datavector, with
		# a spectrum type and a tuple of z bin indices for each C(l)

		# Do a consistency check with the theory vector to ensure the binning is the same
		self.l_bins = block[cl_sections[self.spectra[0][0]], 'ell']
		#if len(self.l_bins) != len(self.data_x):
		#	print 'Error: the angular frequency binning in the theory vector is inconsistent with the data.'
		#	exit()
		if len(self.spectra_theory) != len(self.spectra):
			print 'Error: the selected redshift bin combinations in the theory vector are inconsistent with the data.'
			exit()

		# Set the datavector names now the binning is known
		self.y_name = ['bin_%d_%d'%sp[1] for sp in self.spectra ]

		# Finally evaluate the covariance matrix and invert it
		if not self.constant_covariance:
			self.cov = self.build_covariance(block)
			self.inv_cov = self.build_inverse_covariance(block)

	def build_covariance(self, block):
		"""
		Fits a series of sub-matrices together to give the covariance
		matrix for the full datavector. 
		"""
		cov=None
		if self.constant_covariance:
			cov = self.cov
		else:	
			print 'Calculating covariance matrix.'
			# Get some parameters for this galaxy catalogue needed for the calculation
			self.A_shear = block[self.shear_cat, 'area'] * (np.pi*np.pi)/(180*180)
			self.A_pos = block[self.pos_cat, 'area'] * (np.pi*np.pi)/(180*180)
			self.A = self.A_shear 
			if (self.A_shear!=self.A_pos):
				if self.A_shear<self.A_pos: self.A = self.A_shear
				else: self.A = self.A_pos
			self.l_bins_shear = self.l_bins
			self.l_bins_pos = self.l_bins
			self.dl = block[cl_sections[self.spectra[0][0]], 'l_bin_edges'][1:] - block[cl_sections[self.spectra[0][0]], 'l_bin_edges'][:-1]

		nspec = len(self.spectra)
		nl = self.nlbin
		#n_z_pair = self.nzbin * (self.nzbin+1)/2
		#self.nt = nspec*nl
		nt = self.n_elem
		# Initialise a blank matrix
		Cov = np.zeros((nt, nt))

		# The covariance matrix for the full datavector has a n_lxn_l 
		# block for each pair of spectra

		i,j=0,0
		k = np.zeros(nt)
		# Loop over the types of spectra e.g. ee, nn, ne 
		# This should give a maximum of 9 combinations
		for spect1 in [ x for x in self.spectra if x[1]==(0,0) ]:
			for spect2 in [ y for y in self.spectra if y[1]==(0,0) ]:
				
				m = self.compute_covariance_matrix( block, mode=[ spect1[0],spect2[0]] )
				print 'Evaluated matrix %s%s' %(spect1[0],spect2[0]) 

				# Write to the appropriate block of the matrix
				Cov[ k[i] : (k[i]+m.shape[0]), j : (j+m.shape[1]) ] += m

				# Update indices
				j+=m.shape[1]
				k[i] += m.shape[0]
				i+=1
			i=0
			j=0

		self.cov = Cov
		return Cov

	def normalise_likelihood(self, block):
		sign, logdet = np.linalg.slogdet(self.cov)
		p = 0.5 * logdet
		
		if self.constant_covariance:
			p = 0


		# Read the likelihood, add the normalisation factor and
		# resave to the datablock
		like = block["likelihoods", self.like_name+"_LIKE"]
		like += p
		block.replace_double("likelihoods", self.like_name+"_LIKE", like)

	def build_inverse_covariance(self, block):
		try:
			cov = self.cov
		except:
			self.build_covariance(block)

		self.inv_cov = np.linalg.inv(self.cov)
		block['data_vector', self.like_name+'_inverse_covariance']=self.inv_cov
		return self.inv_cov

	def extract_covariance(self):
		self.cov = np.loadtxt(self.covmat_file)
		print 'Loading covariance matrix from %s' %self.covmat_file

	def build_data(self, datafile):
		print 'Loading datavector.'
		cl_data_all = pickle.load(open(datafile, 'rb'))

		self.data_x = cl_data_all['ell_shear']

		self.nlbin = len(self.data_x)

		self.data_y = []
		self.spectra = []
		for spect in [sp for sp in cl_data_all.keys() if 'C_' in sp]:
			if spect.split('_')[1] not in self.spectra_to_use: continue
			for bin in [b for b in cl_data_all[spect].keys()]:
				b1 = int(bin.split('_')[1]) - 1
				b2 = int(bin.split('_')[2]) - 1
				zcomb = (b1,b2)
				# if the spectrum is symmetric, only load half of the possible bin
				# combinations
				if (spect.split('_')[1][0]==spect.split('_')[1][1]) and (zcomb[1]<zcomb[0]): continue
				else:
					self.spectra.append((spect.split('_')[1], zcomb))
					self.data_y.append( cl_data_all[spect][bin] )
		self.data_y = np.array(self.data_y).flatten()

		return self.data_x, self.data_y

	def extract_theory_points(self, block):
		"""
		Reads the theory spectra from the datablock in the 
		appropriate order. 
		"""
		self.theory_vector = []
		for spect in self.spectra:
			section = cl_sections[spect[0]]
			bin = 'bin_%d_%d' %(spect[1][0]+1, spect[1][1]+1)
			if self.interpolate:
				ell = self.data_x
				theory_vector = np.interp(ell, block[section,'ell'], block[section, bin])
			else:
				ell = block[section, 'ell']
				theory_vector = block[section, bin]
			if self.cuts:
				lmin= block[section, 'lmin_%d_%d'%(spect[1][0]+1, spect[1][1]+1)]
				lmax= block[section, 'lmax_%d_%d'%(spect[1][0]+1, spect[1][1]+1)]
			else:
				print 'No prescription for scale cuts found. Using full datavector.'
				lmin = 0. ; lmax = 1e7

			scale_cut = (ell>lmin) & (ell<lmax)
			theory_vector = theory_vector[scale_cut]

			self.theory_vector.append(theory_vector)

		block['data_vector', self.like_name+'_theory'] = np.concatenate((np.array(self.theory_vector)))
		return np.concatenate((np.array(self.theory_vector)))

	def apply_scale_cuts(self, block, cuts = True):

		self.cuts = cuts
		# Reinitialise the datavector without cuts
		self.data_y = np.copy(self.data_y_uncut)
		cut_data = []
		for i in xrange(len(self.spectra)):
			spect = self.spectra[i]
			section = cl_sections[spect[0]]
			ell = block[section, 'ell']
			if self.cuts:
				lmin= block[section, 'lmin_%d_%d'%(spect[1][0]+1, spect[1][1]+1)]
				lmax= block[section, 'lmax_%d_%d'%(spect[1][0]+1, spect[1][1]+1)]
			else:
				print 'No prescription for scale cuts found. Using the full datavector.'
				return 1
			scale_cut = (ell>lmin) & (ell<lmax)
			
			try: cut_data.append(self.data_y[i*self.nlbin : (i+1)*self.nlbin][scale_cut])
			except: pdb.set_trace()

		self.data_y = np.concatenate((np.array(cut_data)))

		self.nlbin = len(self.data_x)

	def get_zbin_combinations(self, block, auto, cross):
		"""
		Looks up the number of redshift bins and saves all 
		possible combinations of pairs as a list of tuples. 
		"""
		self.bins_ee = []
		self.bins_ne = []
		self.bins_nn = []
		if self.shear:
			self.nzbin_shear = int(block[self.shear_cat, 'nzbin'])

			for i in xrange(self.nzbin_shear):
				for j in xrange(self.nzbin_shear):
					self.bins_ee.append((i,j))

		if self.pos:
			self.nzbin_pos = int(block[self.pos_cat, 'nzbin'])
			for i in xrange(self.nzbin_pos):
				for j in xrange(self.nzbin_pos):
					self.bins_nn.append((i,j))

		if self.ggl:
			self.nzbin_pos = int(block[self.pos_cat, 'nzbin'])
			self.nzbin_shear = int(block[self.shear_cat, 'nzbin'])
			
			for i in xrange(self.nzbin_pos):
				for j in xrange(self.nzbin_shear):
					self.bins_ne.append((i,j))

	def compute_covariance_matrix(self, block, mode):
		""" 
		Calculates an analytic Gaussian covariance matrix for a  
		given combination of spectra and redshift bins using the 
		theory Cls. See e.g. Joachimi and Bridle (2010) eq. (35). 
		"""

		# The covariance involves sevaral cross terms so four type
		# indices need to be extracted, two from each spectrum.
		# mode is a list of two strings e.g. ['ee','ee']
		indices = [(mode[0][0], mode[1][0]), 
			   (mode[0][1], mode[1][1]), 
			   (mode[0][0], mode[1][1]), 
			   (mode[0][1], mode[1][0])]

		# The idea here is to look up the correct sections for each
		# index combination using a dictionary containing section 
		# names
		# The indices may need to be reversed in some cases as the ne 
		# and en spectra are not saved separately, since C^ij_ne = C^ji_en 
		sections=[]
		for i in indices:
			try:
				sections.append( (cl_sections[ '%s%s' %(i[0],i[1]) ], False) )
			except:
				sections.append( (cl_sections[ '%s%s' %(i[1],i[0]) ], True) )

		p = 2. * np.pi / ( self.l_bins * self.dl * self.A ) 

		#import pdb ; pdb.set_trace()
		#n_z_pair = self.nzbin * (self.nzbin+1)/2
		n_z_pair1 = len([i[1] for i in self.spectra if i[0]==mode[0]])
		n_z_pair2 = len([i[1] for i in self.spectra if i[0]==mode[1]])
		cov = np.zeros((n_z_pair1* self.nlbin, n_z_pair2* self.nlbin))

		# Get the l coordinates along each axis of this block, accounting for differences in 
		# redshift binning and scale cuts 
		ell_0=[] 
		ell_1=[]
		ell = block[cl_sections[mode[0]], 'ell']
		for i in self.spectra:
			if i[0]==mode[0]:
				if self.cuts:
					l = np.argwhere((ell > block[cl_sections[mode[0]], 'lmin_%d_%d'%(i[1][0]+1, i[1][1]+1)]) & (ell < block[cl_sections[mode[0]], 'lmax_%d_%d'%(i[1][0]+1, i[1][1]+1)]) ).flatten()
				else:
					l = np.argwhere((ell > 0) & (ell < 1e6 )).flatten()

				for m in l: ell_0 += [(m, i[1])]
		#pdb.set_trace()

		ell = block[cl_sections[mode[1]], 'ell']
		for i in self.spectra:
			if i[0]==mode[1]:
				if self.cuts:
					l = np.argwhere((ell > block[cl_sections[mode[0]], 'lmin_%d_%d'%(i[1][0]+1, i[1][1]+1)]) & (ell < block[cl_sections[mode[0]], 'lmax_%d_%d'%(i[1][0]+1, i[1][1]+1)]) ).flatten()
				else:
					l = np.argwhere((ell > 0) & (ell < 1e6 )).flatten()
				for m in l: ell_1 += [(m, i[1])]

		cov = np.zeros((len(ell_0), len(ell_1)))

		if mode[0]==mode[1]:
			symmetric=True
		else:
			symmetric=False

		out = []
		spect00 = (0,0)
		spect10 = (0,0)
		dat=np.array([])
		#Then just cycle through each point in the flat matrix block
		for m0 in enumerate(ell_0):
			for m1 in enumerate(ell_1):
				y, (l0, spect0) = m0
				x, (l1, spect1) = m1 
				
				# Different l points are trivially 0
				# Also only evaluate half the block in the cases where it's
				# known to be symmetric
				if l0!=l1 or (x>y and symmetric): continue

				# Choose the needed combinations of the four redshift bins
				a0,b0 = spect0[0],spect1[0] # i,k
				a1,b1 = spect0[1],spect1[1] # j,m
				a2,b2 = spect0[0],spect1[1] # i,m
				a3,b3 = spect0[1],spect1[0] # j,k

				for sec in xrange(len(sections)):
					exec 'a,b = a%d,b%d' %(sec,sec) 
					if sections[sec][1]: a,b = b,a
					
					# There's no need to apply scale cuts here as the structure 
					# of ell_0 and ell_1 already take them into account
					exec "Cl%d = block[sections[%d][0], 'bin_%d_%d'][%d]" %(sec,sec,1+a,1+b, l0)
				# Evaluate this point
				m = p[l0] * (Cl0 * Cl1 + Cl2 * Cl3)
				if not np.isfinite(m): pdb.set_trace()

				#pdb.set_trace()

				cov[y][x]+=m

				# If possible evaluate the transpose element by symmetry
				if symmetric:
					if cov[x][y]==0.0:
						cov[x][y]+=m
				spect00 = spect0
				spect10 = spect1

		## For code comparison
		## output data and covariance in a form readable by cl_like
		dat=np.array([])
		data=np.array([])
		for i in range(6):
			sp = np.array([int(self.spectra[i][1][1]+1), int(self.spectra[i][1][0]+1)])
			out = self.data_y[i*12:(i+1)*12]
			out = np.concatenate((sp,out))
			if i==0: data = out
			else: data = np.vstack((data,out))
		#np.savetxt('data.txt', data) 
			

		for i in range(6):
			for j in range(6):
				
				#pdb.set_trace()
				sp = np.array([int(self.spectra[i][1][1]+1), int(self.spectra[i][1][0]+1), int(self.spectra[j][1][1]+1), int(self.spectra[j][1][0]+1)])
				mat = np.diag(cov[j*12:(j+1)*12, i*12:(i+1)*12] )
				mat = np.concatenate((sp, mat))
				if i==0 and j==0: 
					dat=mat
				else:
					dat=np.vstack((dat,  mat ))
		#np.savetxt('covmat.txt', dat)


		return cov




#		n1,n2 = 0,0
#		x0=[0]
#		y0=0
		# Loop over the redshift bin combinations in the right order
#		for spect1 in [i[1] for i in self.spectra if i[0]==mode[0]]:
#			for spect2 in [j[1] for j in self.spectra if j[0]==mode[1]]:
				
				# Only do the calculation for half of the block if possible
				# Each nlxnl sub-section can be filled in above and below
				# the diagonal only if this block is for like correlations
				# e.g. ee,ee ; nn,nn; ne,ne
				
				# Choose the needed combinations of the four redshift bins
#				a0,b0 = spect1[0],spect2[0] # i,k
#				a1,b1 = spect1[1],spect2[1] # j,m
#				a2,b2 = spect1[0],spect2[1] # i,m
#				a3,b3 = spect1[1],spect2[0] # j,k
#				# Evaluate each of the four different spectra, 
				# switching the bin ordering if necessary
#				lmin = block[cl_sections[mode[0]], 'lmin_%d_%d'%(1+spect2[0],1+spect2[1])]
#				lmax = block[cl_sections[mode[0]], 'lmax_%d_%d'%(1+spect2[0],1+spect2[1])]
#				ell = block[cl_sections[mode[0]], 'ell']
#				scale_cut = (ell>lmin) & (ell<lmax)
#				if n2<n1 and (mode[0]==mode[1]): 
#					y0 += len(ell[scale_cut])
#					if n1==0: x0+=[len(ell[scale_cut])]
#					else: x0[n2]+=len(ell[scale_cut])
#					n2+=1
#					continue
#
#				for sec in range(len(sections)):
#					exec 'a,b = a%d,b%d' %(sec,sec) 
#					if sections[sec][1]: a,b = b,a
					
					# Apply THE SAME scale cut to all spectra used to get the covariance for 
					# this particular mode and redshift bin combination
#					exec "Cl%d = block[sections[%d][0], 'bin_%d_%d'][scale_cut]" %(sec,sec,1+a,1+b)
				# Finally allocate the diagonal matrix to the right block of this 
				# covariance matrix
				#pdb.set_trace()
#				m = np.diag( p[scale_cut] * (Cl0 * Cl1 + Cl2 * Cl3) )
				

				# Write to the appropriate part of the matrix for this pair of 
				# spectrum types
#				y,x = m.shape
#				if n1==0:
##					x0+=[x]
#					cov[y0:y0+y, :x] += m
#					trans_cov = cov[:y, y0:y0+x]
#					print 'a writing to: %d:%d, 0:%d' %(y0,y0+y,x) 
#				else:
#					cov[y0:y0+y, x0[n2]:x0[n2]+x] += m
#					trans_cov = cov[x0[n2]:x0[n2]+y, y0:y0+x ]
#					print 'b writing to: %d:%d, %d:%d' %(y0,y0+y,x0[n2],x0[n2]+x)
#
#					np.array(x0)[n2]+=x
#				# If the transpose hasn't been written and this block is symmetric
#				# evaluate the transpose by symmetry
#				if (sum(trans_cov.flatten())==0) and (mode[0]==mode[1]):
#					if n1==0: cov[:y, y0:y0+x]+=m
#					else: cov[x0[n2]:x0[n2]+y, y0:y0+x ]+=m
#				
#				nlbin = len(ell[scale_cut])
#				#pdb.set_trace()
#				print n1, n2 
#
#			#	try:	
#			#		cov[n1*nlbin:(n1+1)*nlbin, n2*nlbin:(n2+1)*nlbin] += m
#			#		trans_cov = cov[n2*nlbin:(n2+1)*nlbin, n1*nlbin:(n1+1)*nlbin ]
#					# If the transpose hasn't been written and this block is symmetric
#			#		if (sum(trans_cov.flatten())==0) and (mode[0]==mode[1]):
#			#			cov[n2*nlbin:(n2+1)*nlbin, n1*nlbin:(n1+1)*nlbin ]+=m
#			#	except: 
#			#		print 'Warning: invalid matrix dimensions'
#			#		pdb.set_trace()
#				#pdb.set_trace()
#				
#				if n1!=0: x0[n2]+=x
#				n2+=1
#				y0+=y
#				import pylab as plt
#				plt.matshow(np.log(cov[:30,:30]))
#				pdb.set_trace()
#			n2 = 0
#			y0=0
#			n1+=1
#
#		return cov
