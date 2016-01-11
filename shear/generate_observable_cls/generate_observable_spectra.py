import scipy.interpolate
import numpy as np

class Cl_class:
	def __init__(self, block, config):

		# Spectra to use
		self.shear = config['shear']
		self.intrinsic_alignments= config['intrinsic_alignments']
		self.clustering = config['clustering']
		self.magnification = config['magnification']

		self.noise = config['noise']
		self.bias = config['bias'][0]
		self.m_per_bin = config['bias'][1]

		shear_cat = config['shear_cat']
		pos_cat = config['pos_cat']

		self.dobinning = config['binning']
		self.window = config['window']

		# Relevant parameters for the noise
		if self.noise: 
			self.sigma_gamma = block.get_double(shear_cat, 'shape_dispersion')
			self.ngal_shear = block.get_double(shear_cat, 'ngal')
			self.ngal_pos = block.get_double(pos_cat, 'ngal')

		# And the configuration of the theory spectra
		self.Nzbin_shear= int(block[shear_cat,'nzbin'])
		self.Nzbin_pos= int(block[pos_cat,'nzbin'])
		self.zbin_edges_shear = [ block[shear_cat, 'edge_%d'%i] for i in xrange(1, self.Nzbin_shear + 1) ]
		self.zbin_edges_pos = [ block[pos_cat, 'edge_%d'%i] for i in xrange(1, self.Nzbin_pos + 1) ]
		if self.shear: 
			if self.intrinsic_alignments:
				gg_section = 'shear_cl_gg'
			else:
				gg_section = 'shear_cl'
			self.l_shear = block[gg_section,'ell']
			self.Nl_shear = len(self.l_shear)
			if self.dobinning: 
				self.Nlbin_shear = int(config['nlbin_shear']) 
			else: 
				self.Nlbin_shear = self.Nl_shear
		if self.clustering: 
			self.l_pos = block['matter_cl','ell']
			self.Nl_pos = len(self.l_pos)
			if self.dobinning: 
				self.Nlbin_pos = int(config['nlbin_pos']) 
			else: 
				self.Nlbin_pos = self.Nl_pos
		if self.shear and self.clustering: 
			self.l_pos = block['matter_cl','ell']
			self.Nl_pos = len(self.l_pos)
			if self.dobinning: 
				self.Nlbin_ggl = int(config['nlbin_ggl']) 
			else: 
				self.Nlbin_ggl = self.Nl_pos

		if self.bias:
			self.multiplicative_bias = [ block[shear_cat, "m%d"%i] for i in xrange(1, self.Nzbin_shear+1) ] 
			self.additive_bias = [ block[shear_cat, "c%d"%i] for i in xrange(1, self.Nzbin_shear+1) ]

		# Finally get the desired binning
		self.get_l_bins(config)

	def load_and_generate_observable_cls(self, block, names):

		# Set up somewhere to put the observable spectra
		if self.shear:
			self.C_ee = np.zeros((self.Nzbin_shear, self.Nzbin_shear, self.Nl_shear))
			self.C_ee_binned = np.zeros((self.Nzbin_shear, self.Nzbin_shear, self.Nlbin_shear))
		if self.clustering:
			self.C_nn = np.zeros((self.Nzbin_pos, self.Nzbin_pos, self.Nl_pos))
			self.C_nn_binned = np.zeros((self.Nzbin_pos, self.Nzbin_pos, self.Nlbin_pos))
		if self.shear and self.clustering:
			self.C_ne = np.zeros((self.Nzbin_pos, self.Nzbin_shear, self.Nl_pos))
			self.C_ne_binned = np.zeros((self.Nzbin_pos, self.Nzbin_shear, self.Nlbin_ggl))
	
		# Then cycle through all the redshift bin combinations
		for i in xrange(1, self.Nzbin_shear+1):
			for j in xrange(1, self.Nzbin_shear+1):
				bin = "bin_%d_%d" %(i,j)
				bin_tr = "bin_%d_%d" %(j,i)

				# The C_GG,II,mm,gg spectra are symmetric
				# This is just bookkeeping to account for the fact we only have half of them
				if (j<i):	a = bin
				else:		a = bin_tr
					
				if self.shear:
											# GG
					if self.intrinsic_alignments:
						self.C_ee[i-1][j-1] += block.get_double_array_1d(names.shear_cl_gg, a)
						self.C_ee[i-1][j-1] += block.get_double_array_1d(names.shear_cl_gi, bin)				# GI
						self.C_ee[i-1][j-1] += block.get_double_array_1d(names.shear_cl_gi, bin_tr)					# IG
						try:	self.C_ee[i-1][j-1] += block.get_double_array_1d(names.shear_cl_ii, a)					# I
						except: print "WARNING: No II spectrum found."
					else:
						self.C_ee[i-1][j-1] += block.get_double_array_1d(names.shear_cl, a)

		for i in xrange(1, self.Nzbin_pos+1):
			for j in xrange(1, self.Nzbin_pos+1):
				bin = "bin_%d_%d" %(i,j)
				bin_tr = "bin_%d_%d" %(j,i)

				# The C_GG,II,mm,gg spectra are symmetric
				# This is just bookkeeping to account for the fact we only have half of them
				if (j<i):	a = bin
				else:		a = bin_tr

				if self.clustering:			
					self.C_nn[i-1][j-1] += block.get_double_array_1d('matter_cl', a)							# gg
					if self.magnification:
						self.C_nn[i-1][j-1] += block.get_double_array_1d(names.galaxy_magnification_cl, bin)				# mg
						self.C_nn[i-1][j-1] += block.get_double_array_1d(names.galaxy_magnification_cl, bin_tr)			# gm
						self.C_nn[i-1][j-1] += block.get_double_array_1d(names.magnification_magnification_cl, a)			# mm
		for i in xrange(1, self.Nzbin_pos+1):
			for j in xrange(1, self.Nzbin_shear+1):
				bin = "bin_%d_%d" %(i,j)
				bin_tr = "bin_%d_%d" %(j,i)

				# The C_GG,II,mm,gg spectra are symmetric
				# This is just bookkeeping to account for the fact we only have half of them
				if (j<i):	a = bin
				else:		a = bin_tr
				if self.shear and self.clustering:
					try: 
						self.C_ne[i-1][j-1] += block.get_double_array_1d(names.ggl_cl, bin)							# gG
						if self.intrinsic_alignments:
							self.C_ne[i-1][j-1] += block.get_double_array_1d(names.gal_IA_cross_cl, bin)					# gI
							if self.magnification:
								self.C_ne[i-1][j-1] += block.get_double_array_1d(names.magnification_intrinsic_cl, bin)		# mI
						if self.magnification:
							self.C_ne[i-1][j-1] += block.get_double_array_1d(names.magnification_shear_cl, bin)				# mG
					except:  print "WARNING: One of the components of the GGL spectrum is missing."
	
				if not self.noise:
					# Finally resample the spectra in the survey angular frequency bins
					if self.shear:
						self.C_ee_binned[i-1][j-1] = get_binned_cl(self.C_ee[i-1][j-1], self.l_shear, self.lbin_edges_shear, self.dobinning, self.window )
					if self.clustering:
						self.C_nn_binned[i-1][j-1] = get_binned_cl(self.C_nn[i-1][j-1], self.l_pos, self.lbin_edges_pos, self.dobinning, self.window )
					if self.shear and self.clustering:
						self.C_ne_binned[i-1][j-1] = get_binned_cl(self.C_ne[i-1][j-1], self.l_pos, self.lbin_edges_pos, self.dobinning, self.window )

		if self.noise:
			# Add shot noise if required	
			self.add_noise(block)
			# If noise was added earlier, the binning is done here rather than 
			# immediately on loading
			if self.shear:
				for i in xrange(1, self.Nzbin_shear+1):
					for j in xrange(1, self.Nzbin_shear+1):
						self.C_ee_binned[i-1][j-1] = get_binned_cl(self.C_ee[i-1][j-1], self.l_shear, self.lbin_edges_shear, self.dobinning, self.window )
						if self.bias: self.apply_measurement_bias(i, j, 'shear')
			if self.clustering:
				for i in xrange(1, self.Nzbin_pos+1):
					for j in xrange(1, self.Nzbin_pos+1):
						self.C_nn_binned[i-1][j-1] = get_binned_cl(self.C_nn[i-1][j-1], self.l_pos, self.lbin_edges_pos, self.dobinning, self.window )
			if self.shear and self.clustering:
				for i in xrange(1, self.Nzbin_pos+1):
					for j in xrange(1, self.Nzbin_shear+1):
						self.C_ne_binned[i-1][j-1] = get_binned_cl(self.C_ne[i-1][j-1], self.l_pos, self.lbin_edges_pos, self.dobinning, self.window )
						if self.bias: self.apply_measurement_bias(i, j,'ggl')

	def apply_measurement_bias(self, i, j, mode=None):
		if not self.m_per_bin:
			m0 = self.multiplicative_bias[0]

		#Compute scaling parameter for this pair of redshift bins
		if self.m_per_bin:
			mi = self.multiplicative_bias[i-1]
			mj = self.multiplicative_bias[j-1]
		else:
			mi,mj = m0,m0

		#Apply scaling
		if mode=='shear': 
			self.C_ee_binned[i-1][j-1] *= (1+mi)*(1+mj)
		if mode=='ggl': 
			self.C_ne_binned[i-1][j-1] *= (1+mj)
		
	def add_noise(self, block):
		
		n_binned_shear = get_binned_number_densities(self.Nzbin_shear, self.ngal_shear)
		n_binned_pos = get_binned_number_densities(self.Nzbin_pos, self.ngal_pos)

		# Create noise matrices with the same shape as the Cls
		# These are diagonal in the x,z plane (fixed l) and constant along the y axis (constant redshift)

		N_ee_0 = np.identity(self.Nzbin_shear) * self.sigma_gamma**2 / (2. * n_binned_shear)
		N_nn_0 = np.identity(self.Nzbin_pos) * 1. / n_binned_pos

		N_shot_ee = [] ; N_shot_nn = []

		for i in xrange( len(self.C_ee[0][0]) ):
			N_shot_ee += [ N_ee_0 ]
			N_shot_nn += [ N_nn_0 ]

		N_shot_nn = np.swapaxes(N_shot_nn,0,2) ; N_shot_nn = np.swapaxes(N_shot_nn,0,1) 
		N_shot_ee = np.swapaxes(N_shot_ee,0,2) ; N_shape = np.swapaxes(N_shot_ee,0,1)

		# Then add the relevant noise to the Cl matrices
		if self.shear:	self.C_ee += N_shot_ee
		if self.clustering:	self.C_nn += N_shot_nn

	def get_l_bins(self, config):

		if self.dobinning:
			# Define some l bins for these galaxy samples
			lmin, lmax= config['lmin_shear'], config['lmax_shear']
			if self.shear:
				self.lbin_edges_shear = np.logspace(np.log10(lmin), np.log10(lmax), self.Nlbin_shear+1)
				self.l_bins_shear = np.zeros_like(self.lbin_edges_shear[:-1])
				if self.window == 'tophat':
					self.l_bins_shear = np.exp( (np.log(self.lbin_edges_shear[1:] * self.lbin_edges_shear[:-1]))/2.0 ) 
				elif self.window == 'delta':
					for i in xrange(len(self.lbin_edges_shear)-1):
						lmin0 = self.lbin_edges_shear[i]
						lmax0 = self.lbin_edges_shear[i+1]
						sel = (self.l_shear>lmin0) & (self.l_shear<lmax0)
						l_in_window = self.l_shear[sel]
						self.l_bins_shear[i] = l_in_window[len(l_in_window)/2]

			if self.clustering:
				self.lbin_edges_pos = np.logspace(np.log10(lmin), np.log10(lmax), self.Nlbin_pos+1)
				self.l_bins_pos = np.zeros_like(self.lbin_edges_pos[:-1])
				lmin, lmax= config['lmin_pos'], config['lmax_pos']
				if self.window == 'tophat':
					self.l_bins_pos = np.exp( (np.log(self.lbin_edges_pos[1:] * self.lbin_edges_pos[:-1]))/2.0 ) 
				elif self.window == 'delta':
					for i in xrange(len(self.lbin_edges_pos)-1):
						lmin0 = self.lbin_edges_pos[i]
						lmax0 = self.lbin_edges_pos[i+1]
						sel = (self.l_pos>lmin0) & (self.l_pos<lmax0)
						l_in_window = self.l_pos[sel]
						self.l_bins_pos[i] = l_in_window[len(l_in_window)/2]
			
			if self.clustering and self.shear:
				self.lbin_edges_ggl = np.logspace(np.log10(lmin), np.log10(lmax), self.Nlbin_ggl+1)
				self.l_bins_ggl = np.zeros_like(self.lbin_edges_ggl[:-1])
				lmin, lmax= config['lmin_ggl'], config['lmax_ggl']
				self.lbin_edges_ggl = np.logspace(np.log10(lmin), np.log10(lmax), self.Nlbin_ggl+1)
				if self.window == 'tophat':
					self.l_bins_ggl = np.exp( (np.log(self.lbin_edges_ggl[1:] * self.lbin_edges_ggl[:-1]))/2.0 ) 
				elif self.window == 'delta':
					for i in xrange(len(self.lbin_edges_ggl)-1):
						lmin0 = self.lbin_edges_ggl[i]
						lmax0 = self.lbin_edges_ggl[i+1]
						sel = (self.l_ggl>lmin0) & (self.l_ggl<lmax0)
						l_in_window = self.l_ggl[sel]
						self.l_bins_ggl[i] = l_in_window[len(l_in_window)/2]
		else:
			if self.shear:
				self.l_bins_shear = self.l_shear
				self.lbin_edges_shear = None
			if self.clustering:
				self.l_bins_pos = self.l_pos
				self.lbin_edges_pos = None
			if self.shear and self.clustering:
				self.l_bins_ggl = self.l_ggl
				self.lbin_edges_ggl = None
			


	def save_cls(self, block, out_path):
		
		if self.shear:
			if self.dobinning:
				block.put_double_array_1d('galaxy_shape_cl', 'l_bin_edges', self.lbin_edges_shear)
			block.put_double_array_1d('galaxy_shape_cl', 'ell', self.l_bins_shear)
			block.put_double_array_1d('galaxy_shape_cl', 'z_bin_edges', self.zbin_edges_shear)
			block.put_int('galaxy_shape_cl', 'nl', self.Nlbin_shear)
			block.put_int('galaxy_shape_cl', 'nz', self.Nzbin_shear)
		if self.clustering:
			if self.dobinning:
				block.put_double_array_1d('galaxy_position_cl', 'l_bin_edges', self.lbin_edges_pos)
			block.put_double_array_1d('galaxy_position_cl', 'ell', self.l_bins_pos)
			block.put_double_array_1d('galaxy_position_cl', 'z_bin_edges', self.zbin_edges_pos)
			block.put_int('galaxy_position_cl', 'nl', self.Nlbin_pos)
			block.put_int('galaxy_position_cl', 'nz', self.Nzbin_pos)
		if self.shear and self.clustering:
			if self.dobinning:
				block.put_double_array_1d('galaxy_position_shape_cross_cl', 'l_bin_edges', self.lbin_edges_ggl)
			block.put_double_array_1d('galaxy_position_shape_cross_cl', 'ell', self.l_bins_ggl)
			block.put_double_array_1d('galaxy_position_shape_cross_cl', 'z_bin_edges_shear', self.zbin_edges_shear)
			block.put_double_array_1d('galaxy_position_shape_cross_cl', 'z_bin_edges_position', self.zbin_edges_pos)
			block.put_int('galaxy_position_shape_cross_cl', 'nl', self.Nlbin_ggl)
			block.put_int('galaxy_position_shape_cross_cl', 'nz_shear', self.Nzbin_shear)			
			block.put_int('galaxy_position_shape_cross_cl', 'nz_position', self.Nzbin_pos)	

		out_C_ee={}
		out_C_nn={}
		out_C_ne={}
		if self.shear:
			for i in xrange(1, self.Nzbin_shear+1):
				for j in xrange(1, self.Nzbin_shear+1):
					bin = "bin_%d_%d" %(i,j)
					block.put_double_array_1d('galaxy_shape_cl', bin, self.C_ee_binned[i-1][j-1])
					if out_path is not None:	out_C_ee[bin]=	self.C_ee_binned[i-1][j-1]
		if self.clustering:
			for i in xrange(1, self.Nzbin_pos+1):
				for j in xrange(1, self.Nzbin_pos+1):
					bin = "bin_%d_%d" %(i,j)
					block.put_double_array_1d('galaxy_position_cl', bin, self.C_nn_binned[i-1][j-1])
					if out_path is not None:	out_C_nn[bin]=	self.C_nn_binned[i-1][j-1]
		if self.shear and self.clustering:
			for i in xrange(1, self.Nzbin_pos+1):
				for j in xrange(1, self.Nzbin_shear+1):
					bin = "bin_%d_%d" %(i,j)
					block.put_double_array_1d('galaxy_position_shape_cross_cl', bin, self.C_ne_binned[i-1][j-1])
					if out_path is not None:	out_C_ne[bin]=	self.C_ne_binned[i-1][j-1]
		
		try: omega_de = block['cosmological_parameters', 'omega_de']
		except:
			omega_k = block['cosmological_parameters', 'omega_k']
			omega_de = 1.0 - block['cosmological_parameters', 'omega_m'] - omega_k
		if out_path is not None:
			cospar={'omega_m': block['cosmological_parameters', 'omega_m'],
				'omega_de': omega_de,
				'omega_b': block['cosmological_parameters', 'omega_b'],
				'h': block['cosmological_parameters', 'h0'],
				'sigma_8': block['cosmological_parameters', 'sigma_8'],
				'n_s': block['cosmological_parameters', 'n_s'],
				'w0': block['cosmological_parameters', 'w'],
				'wa': block['cosmological_parameters', 'wa']
				}
			datavector = {'cosmology': cospar}

			if self.shear: 
				datavector['C_ee'] = out_C_ee
				datavector['l_bin_edges_shear'] = self.lbin_edges_shear
				datavector['z_bin_edges_shear'] = self.zbin_edges_shear
				datavector['ell_shear'] = self.l_bins_shear
			if self.clustering: 
				datavector['C_nn'] = out_C_nn
				datavector['l_bin_edges_pos'] = self.lbin_edges_pos
				datavector['z_bin_edges_pos'] = self.zbin_edges_pos
				datavector['ell_pos'] = self.l_bins_pos
			if self.shear and self.clustering: 
				datavector['C_ne'] = out_C_ne
				datavector['l_bin_edges_ggl'] = self.lbin_edges_ggl
				datavector['ell_ggl'] = self.l_bins_ggl
			import cPickle as pickle

			f = open(out_path, 'wb')
			print 'Saved datavector to %s'%out_path
			pickle.dump(datavector, f)
			f.close()

def get_binned_cl(Cl, l, lbin_edges, dobinning, window):
	if dobinning:	
		Cl_binned = np.zeros(len(lbin_edges)-1)
		for i in xrange(len(lbin_edges)-1):
			lmin = lbin_edges[i] ; lmax = lbin_edges[i+1]
			sel = (l>lmin) & (l<lmax)
			if window == 'tophat':
				Cl_binned[i] = np.mean(Cl[sel])
			elif window == 'delta':
				i0 = len(Cl[sel])/2
				Cl_binned[i] = Cl[sel][i0]
		return Cl_binned
	else:
		return Cl


def get_binned_number_densities(nzbin, ngal):
	"""
	Calculate the average number density of galaxies in each redshift bin,
	assuming an equal number in each. 
	"""
	n_binned = []
	# Convert number density from arcmin^-2 to sr^-1
	ngal = ( 60*60 * 180*180 / (np.pi*np.pi) ) * ngal 		
	for i in xrange(nzbin):
		n_binned += [ngal/nzbin]
	n_binned = np.array(n_binned)

	return n_binned


		
