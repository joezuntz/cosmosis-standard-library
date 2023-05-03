from numpy import log, pi, interp, where, loadtxt, dot, append, linalg, genfromtxt
import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.gaussian_likelihood import GaussianLikelihood

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
growthparams = section_names.growth_parameters
dist = section_names.distances

ROOT_DIR = os.path.split(os.path.abspath(__file__))[0]

c_km_per_s = 299792.458
default_rd_fiducial = 147.78


class Boss12LRGLikelihood(GaussianLikelihood):
	
	data_type = "LRG"
	like_name = "boss12_lrg"
	def __init__(self, options):
		
		super(Boss12LRGLikelihood, self).__init__(options)
		# Allow override of these parameters
		self.rd_fiducial = self.options.get_double("rd_fiducial", default_rd_fiducial)
		self.mode = self.options.get_int("mode", default=0)
		self.feedback = self.options.get_bool("feedback", default=False)
		
	def build_data(self):
		
		self.mode = self.options.get_int("mode", default=0)
		
		if not self.mode:
			print("LRG data")
			print("BAO only: Dm(z)/rs and Dh(z)/rs")
			# Reading data file
			data_file = os.path.join(ROOT_DIR, "sdss_DR12_LRG_BAO_DMDH.txt")
		else:
			print("LRG data")
			print("BAO+FS: Dm(z)/rd, Dh(z)/rd and f(z)sigma8(z)")
			# Reading data file
			data_file = os.path.join(ROOT_DIR, "sdss_DR12_LRG_FSBAO_DMDHfs8.txt")
			
		DATA = loadtxt(data_file, usecols=(0, 1))
		z_eff, data = DATA[:, 0], DATA[:, 1]
		
		return z_eff, data
		
	def build_covariance(self):
		
		if not self.mode:
			# Reading covariance matrix file
			cov_file = os.path.join(ROOT_DIR, 'sdss_DR12_LRG_BAO_DMDH_covtot.txt')
		else:
			# Reading covariance matrix file
			cov_file = os.path.join(ROOT_DIR, 'sdss_DR12_LRG_FSBAO_DMDHfs8_covtot.txt')
			
		cov = loadtxt(cov_file)
		self.inv_cov = linalg.inv(cov)
		
		return cov
		
	def build_inverse_covariance(self):
		return self.inv_cov
		
	def extract_theory_points(self, block):
		
		# Redshift array
		z = block[dist, 'z']
		# Sound horizon at the drag epoch
		rd = block[dist, "rs_zdrag"]
		# Comoving distance
		Dm = block[dist, 'd_m']  # in Mpc
		# Hubble distance
		H = c_km_per_s*block[dist, 'H']  # in c/Mpc
		Dh = c_km_per_s/H
		
		# z and distance maybe are loaded in chronological order
		# Reverse to start from low z
		if (z[1] < z[0]):
			z = z[::-1]
			Dm = Dm[::-1]
			Dh = Dh[::-1]
			
		#Find theory Dm and Dh at effective redshift by interpolation
		z_eff = self.data_x
		
		# Add fsigma8 when mode == 1 (BAO+FS)
		if self.mode:
			Dm_z_rd = interp(z_eff[0::3], z, Dm)/rd
			Dh_z_rd = interp(z_eff[1::3], z, Dh)/rd
			z = block['growth_parameters', 'z']
			fsigma8 = block['growth_parameters', 'fsigma_8']
			# Find theory fsigma8 at fiducial redshift
			fsigma8_z = interp(z_eff[2::3], z, fsigma8)
			params = [Dm_z_rd[0],Dh_z_rd[0],fsigma8_z[0],Dm_z_rd[1],Dh_z_rd[1],fsigma8_z[1]]

			
		else:
			# Distances over rd
			Dm_z_rd = interp(z_eff[0::2], z, Dm)/rd
			Dh_z_rd= interp(z_eff[1::2], z, Dh)/rd
			params = [Dm_z_rd[0],Dh_z_rd[0],Dm_z_rd[1],Dh_z_rd[1]]
			
		if self.feedback:
			print()
			print('             zeff   pred    data')
			if self.mode:
				print('Dm_over_rd: %.3f  %.3f  %.3f' % (self.data_x[0], Dm_z_rd[0], self.data_y[0]))
				print('Dh_over_rd: %.3f  %.3f  %.3f' %(self.data_x[1], Dh_z_rd[0], self.data_y[1]))
				print('fsigma8:    %.3f   %.3f   %.3f' %(self.data_x[2], fsigma8_z[0], self.data_y[2]))
				print('Dm_over_rd: %.3f  %.3f  %.3f' % (self.data_x[3], Dm_z_rd[1], self.data_y[3]))
				print('Dh_over_rd: %.3f  %.3f  %.3f' %(self.data_x[4], Dh_z_rd[1], self.data_y[4]))	
				print('fsigma8:    %.3f   %.3f   %.3f' %(self.data_x[5], fsigma8_z[1], self.data_y[5]))
			else:
				print('Dm_over_rd: %.3f  %.3f  %.3f' % (self.data_x[0], Dm_z_rd[0], self.data_y[0]))
				print('Dh_over_rd: %.3f  %.3f  %.3f' %(self.data_x[1], Dh_z_rd[0], self.data_y[1]))
				print('Dm_over_rd: %.3f  %.3f  %.3f' % (self.data_x[2], Dm_z_rd[1], self.data_y[2]))
				print('Dh_over_rd: %.3f  %.3f  %.3f' %(self.data_x[3], Dh_z_rd[1], self.data_y[3]))	
			print()
			
		return params
		
setup, execute, cleanup = Boss12LRGLikelihood.build_module()
