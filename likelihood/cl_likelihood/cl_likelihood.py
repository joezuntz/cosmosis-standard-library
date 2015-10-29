from gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
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
	def __init__(self, block, survey, covmat_file, data_file, shear, galaxy_clustering, sh_gal_cross, auto=None, cross=None):

		# Put all of the information specifying the datavector where
		# it can be accessed later
		spectra_to_use = []
		if shear:	spectra_to_use.append('ee')
		if sh_gal_cross:	spectra_to_use.append('ne')
		if galaxy_clustering:	spectra_to_use.append('nn')

		self.survey = survey

		self.llim = (block[self.survey, 'lmin'],block[self.survey, 'lmax'])
		bins = self.get_zbin_combinations(block, auto, cross)

		self.spectra = []
		for sp in spectra_to_use:
			for zcomb in bins: 
				# For the symmetric spectra ee and nn, only half the redshift
				# bin combinations are required
				if (sp[0]==sp[1]) and (zcomb[1]<zcomb[0]): continue
				else: self.spectra.append((sp, zcomb))
		# This list should now fully specify the datavector, with
		# a spectrum type and a tuple of z bin indices for each C(l)

		# Load the l bins for the theory vector now, so they can be verified later
		self.l_bins = block[cl_sections[self.spectra[0][0]], 'l_bins']

		self.x_section = [cl_sections[sp] for sp in spectra_to_use]
		self.x_name    = 'l_bin_edges'
		self.y_section = [cl_sections[sp] for sp in spectra_to_use]
		self.y_name = ['bin_%d_%d'%sp[1] for sp in self.spectra ]
		if len(spectra_to_use)==1:
			names.data_vector = "cl_%s"%spectra_to_use[0]
			self.like_name = "cl_%s"%spectra_to_use[0]
		else:
			names.data_vector = 'combined_cl_%s'%string.join(spectra_to_use).replace(' ','_')
			self.like_name = 'combined_cl_%s'%string.join(spectra_to_use).replace(' ','_')

		self.data_x, self.data_y = self.build_data(data_file)
		self.cov = self.build_covariance(block, covmat_file)
		self.inv_cov = self.build_inverse_covariance(block)

	def build_covariance(self, block, filename=None):
		"""
		Fits a series of sub-matrices together to give the covariance
		matrix for the full datavector. 
		"""
		cov = None
		if (filename!=None):
			cov_file = open(filename, 'rb')
			cov = pickle.load(cov_file)
			cov_file.close()
			print 'Loading covariance matrix from %s' %filename
		else:	
			print 'Calculating covariance matrix.'
			# Get some survey parameters needed for the calculation
			self.A = block[self.survey, 'area'] * (np.pi*np.pi)/(180*180)
			self.l_bins = self.l_bins[self.scale_cut]
			self.dl = block[cl_sections[self.spectra[0][0]], 'l_bin_edges'][1:] - block[cl_sections[self.spectra[0][0]], 'l_bin_edges'][:-1]
			self.dl = self.dl[self.scale_cut]		

		nspec = len(self.spectra)
		nl = self.nlbin
		n_z_pair = self.nzbin * (self.nzbin+1)/2
		self.nt = nspec*nl
		nt = self.nt
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
				if (filename!=None):
					m = cov[ 'cov_%s%s'%(spect1[0],spect2[0]) ][spect1[1]][spect2[1]][self.scale_cut]
				else:
					m = self.compute_covariance_matrix( block, mode=[ spect1[0],spect2[0]] )

					print 'Evaluated matrix %s%s' %(spect1[0],spect2[0]) 

				# Write to the appropriate block of the matrix
				try:	Cov[ k[i] : (k[i]+m.shape[0]), j : (j+m.shape[1]) ] += m
				except:	pdb.set_trace()

				#if (spect1[0],spect2[0])== ('nn','nn'): pdb.set_trace()
				# Update indices
				j+=m.shape[1]
				k[i] += m.shape[0]
				i+=1
			i=0
			j=0
		#import matplotlib.pyplot as plt
		#plt.imshow(Cov[0:50,0:50], interpolation='None',origin='lower')

		self.cov = Cov
		return Cov

	def build_inverse_covariance(self, block):
		try:
			cov = self.cov
		except:
			self.build_covariance(block)
		self.inv_cov = np.linalg.inv(self.cov)
		return self.inv_cov

	def build_data(self, datafile):
		print 'Loading datavector.'
		cl_data_all = pickle.load(open(datafile, 'rb'))

		self.data_x = (cl_data_all['l_bin_edges'][1:] + cl_data_all['l_bin_edges'][:-1])/2.
		
		# Do a consistency check with the theory vector to ensure the binning is the same
		if len(self.l_bins) != len(self.data_x):
			print 'Error: the angular frequency binning in the theory vector is inconsistent with the data.'
			exit()


		if (self.llim[0]!=None) and (self.llim[1]!=None): 
			self.scale_cut = (self.data_x>self.llim[0]) & (self.data_x<self.llim[1])
		else:
			self.scale_cut = np.array([True]*len(self.data_x))

		self.data_x = self.data_x[self.scale_cut]
		self.nlbin = len(self.data_x)

		self.data_y = []

		for spect in self.spectra:
			bin = 'bin_%d_%d' %(spect[1][0]+1, spect[1][1]+1)
			self.data_y.append( cl_data_all['C_%s'%spect[0]][bin][self.scale_cut] )
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
			self.theory_vector.append(block[section, bin][self.scale_cut])

		return np.array(self.theory_vector).flatten()

	def get_zbin_combinations(self, block, auto, cross):
		"""
		Looks up the number of redshift bins and saves all 
		possible combinations of pairs as a list of tuples. 
		"""
		self.nzbin = int(block[self.survey, 'nzbin'])
		bins = []

		for i in range(self.nzbin):
			for j in range(self.nzbin):
				bins.append((i,j))

		return bins

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

		n_z_pair = self.nzbin * (self.nzbin+1)/2
		n_z_pair1 = len([i[1] for i in self.spectra if i[0]==mode[0]])
		n_z_pair2 = len([i[1] for i in self.spectra if i[0]==mode[1]])
		cov = np.zeros((n_z_pair1* self.nlbin, n_z_pair2* self.nlbin))

		n1,n2 = 0,0
		# Loop over the redshift bin combinations in the right order
		for spect1 in [i[1] for i in self.spectra if i[0]==mode[0]]:
			for spect2 in [j[1] for j in self.spectra if j[0]==mode[1]]:
				# If possible get the sub-matrix by symmetry
				if (mode[0]==mode[1]) and (spect1>spect2):
					m = cov[n2*self.nlbin:(n2+1)*self.nlbin, n1*self.nlbin:(n1+1)*self.nlbin ]
				# Otherwise calculate it
				else:
					# Choose the needed combinations of the four redshift bins
					a0,b0 = spect1[0],spect2[0] # i,k
					a1,b1 = spect1[1],spect2[1] # j,m
					a2,b2 = spect1[0],spect2[1] # i,m
					a3,b3 = spect1[1],spect2[0] # j,k
					# Evaluate each of the four different spectra, 
					# switching the bin ordering if necessary
					for sec in range(len(sections)):
						exec 'a,b = a%d,b%d' %(sec,sec) 
						if sections[sec][1]: a,b = b,a
						exec "Cl%d = block[sections[%d][0], 'bin_%d_%d'][self.scale_cut]" %(sec,sec,1+a,1+b)
						if mode[1]=='ne': 
							print "Cl%d = block[%s, 'bin_%d_%d'][self.scale_cut]" %(sec,sections[sec][0],1+a,1+b)
					# Finally allocate the diagonal matrix to the right block of this 
					# covariance matrix
					m = np.diag( p * (Cl0 * Cl1 + Cl2 * Cl3) )
				# Write to the appropriate part of the matrix for this pair of 
				# spectrum types
				try:	cov[n1*self.nlbin:(n1+1)*self.nlbin, n2*self.nlbin:(n2+1)*self.nlbin] += m
				except: 
					print 'Warning: invalid matrix dimensions'
					pdb.set_trace()
				
				n2+=1
			n2 = 0
			n1+=1

		return cov
