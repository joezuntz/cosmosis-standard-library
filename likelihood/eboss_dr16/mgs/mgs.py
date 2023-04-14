from numpy import log, pi, interp, where, loadtxt, dot, append, linalg, genfromtxt
import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.gaussian_likelihood import GaussianLikelihood

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
growthparams = section_names.growth_parameters
dist = section_names.distances

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]

c_km_per_s = 299792.458
default_rd_fiducial = 147.8


class MGSLikelihood(GaussianLikelihood):
	
	data_type = "MGS"
	like_name = "mgs"
	def __init__(self, options):
		
		super(MGSLikelihood, self).__init__(options)
		# Allow override of these parameters
		self.rd_fiducial = self.options.get_double("rd_fiducial", default_rd_fiducial)
		self.feedback = self.options.get_bool("feedback", default=False)
		
	def build_data(self):
		
		print("MGS data")
		print("BAO+FS: Dv(z) and f(z)sigma8(z)")
		# Reading data file
		DATA_file = os.path.join(ROOT_dir, "sdss_MGS_FSBAO_DVfs8.txt")
			
		DATA = loadtxt(DATA_file, usecols=(0, 1))
		z_eff, data = DATA[:, 0], DATA[:, 1]
		
		return z_eff, data
		
	def build_covariance(self):
		
		# Reading covariance matrix file
		COV_file = os.path.join(ROOT_dir, 'sdss_MGS_FSBAO_DVfs8_covtot.txt')
			
		cov = loadtxt(COV_file)
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
    	
		# Spherically averaged distance
		Dv = (z * Dh *Dm**2)**(1/3)
		
		#Find theory Dm and Dh at effective redshift by interpolation
		z_eff = self.data_x
		
		Dv_z_rd = interp(z_eff[1], z, Dv)/rd
		z = block['growth_parameters', 'z']
		fsigma8 = block['growth_parameters', 'fsigma_8']
		# Find theory fsigma8 at fiducial redshift
		fsigma8_z = interp(z_eff[0], z, fsigma8)
		params = [fsigma8_z,Dv_z_rd]
		
		if self.feedback:
			print()
			print('             zeff   pred    data')
			print('fsigma8:    %.3f   %.3f   %.3f' %(self.data_x[0], fsigma8_z, self.data_y[0]))
			print('Dv_over_rd: %.3f  %.3f  %.3f' % (self.data_x[1], Dv_z_rd, self.data_y[1]))
			print()
			
		return params
		
setup, execute, cleanup = MGSLikelihood.build_module()
