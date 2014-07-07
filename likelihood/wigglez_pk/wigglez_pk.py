import os
import numpy as np

REGION_SIZE = 50
P_COLUMN = 3
NREGION = 7
import pylab

class WigglezDataset(object):
	def __init__(self, common_filename, dataset_filename, root):
		common_options = self.parse_dataset_file(common_filename)
		options = self.parse_dataset_file(dataset_filename)

		kmin_index = int(common_options['min_mpk_kbands_use']) - 1
		kmax_index = int(common_options['max_mpk_kbands_use']) - 1

		self.name = options["name"]
		self.z = float(options["redshift"])
		self.dv_fid = float(options['DV_fid'])
		self.first_datapoint = int(common_options['min_mpk_points_use']) - 1
		self.last_datapoint = int(common_options['max_mpk_points_use']) - 1

		#Paths to the files
		windows_file = os.path.join(root, options['windows_file'])
		kbands_file = os.path.join(root, options['kbands_file'])
		measurements_file = os.path.join(root, options['measurements_file'])
		cov_file = os.path.join(root, options['cov_file'])

		#Load in the data from the files
		self.window_matrix = np.loadtxt(windows_file).T
		self.k = np.loadtxt(kbands_file)[kmin_index:kmax_index+1]
		self.data = np.loadtxt(measurements_file).T
		self.cov = np.loadtxt(cov_file).T
		self.weight = []
		for i in xrange(NREGION):
			cov_chunk = self.cov[:, i*REGION_SIZE:(i+1)*REGION_SIZE]
			cov_chunk = cov_chunk[self.first_datapoint:self.last_datapoint+1, self.first_datapoint:self.last_datapoint+1]
			w = np.linalg.inv(cov_chunk)
			self.weight.append(w)

		self.region_size = REGION_SIZE

		print "Loaded %s "% self.name
		# print "    Window matrix shape = ", self.window_matrix.shape
		# print "    k shape = ", self.k.shape
		# print "    data shape = ", self.data.shape
		# print "    cov shape = ", self.cov.shape

		print "Not sure if I have the a_scaling part of this code the right way"
		print "around yet.  Do not use in production until this is figured out!"


	def region_likelihood(self, k_theory, p_theory, region):
		chunk = slice(region*REGION_SIZE,(region+1)*REGION_SIZE)
		k_grid = self.k  #This is the convolution k grid, not the obs one
		p_grid = np.interp(k_grid, k_theory, p_theory) #theory P at k grid points

		window = self.window_matrix[:, chunk]
		p_theory = np.dot(p_grid, window)
		p_obs = self.data[P_COLUMN, chunk]
		k_obs = self.data[0, chunk]
		weight = self.weight[region]
		# covmat = self.cov[:, chunk]
		# sigma2 = covmat.diagonal()
		# sigma = sigma2**0.5
		# pylab.subplot(2,4,region+1)
		# pylab.plot(k_obs, p_theory,'-')
		# pylab.errorbar(k_obs, p_obs, sigma, fmt='.')
		# pylab.gca().set_xscale("log", nonposx='clip')
		# pylab.gca().set_yscale("log", nonposy='clip')
		# if region==6: pylab.show()
		# import sys
		# sys.exit()
		# print p_obs-p_theory
		d = (p_obs - p_theory)[self.first_datapoint:self.last_datapoint+1]
		# print d**2/sigma2, (d**2/sigma**2).sum()
		chi2 = np.dot(d, np.dot(weight, d))
		return -0.5*chi2


	def likelihood(self, H0, z1, dv, k, z2, p_k):


		#Interpolate into dv(z1) at my redshift to get my D_V
		dv_theory = np.interp([self.z], z1, dv)  * H0
		a_scaling = (dv_theory / self.dv_fid)**(1./3.)
		# print "a_scaling = ", a_scaling
		# print "D_V theory = ",  dv_theory

		#Interpolate into the 2D p_k(k, z2) to get the 
		#linear slice through at my redshift
		z_index = np.searchsorted(z2, self.z)
		zlow = z2[z_index]
		zhigh = z2[z_index+1]
		dz = zhigh-zlow
		r = (self.z-zlow)/dz
		p_k = p_k[z_index,:]*(1-r) + p_k[z_index+1,:]*r
		
		L = 0.0
		for i in xrange(NREGION):
			L += self.region_likelihood(k*a_scaling, p_k/a_scaling**3, i)
		return L


	@staticmethod
	def parse_dataset_file(filename):
		output = {}
		for line in open(filename):
			line=line.strip()
			if not line or line.startswith('#'): continue
			key, value = line.split('=')
			output[key.strip()] = value.strip()
		return output

def test():
	dist_dirname='../../../demo_output_1/distances/'
	z1 = np.loadtxt(dist_dirname+"z.txt")
	d_a = np.loadtxt(dist_dirname+"d_a.txt")
	h = np.loadtxt(dist_dirname+"h.txt")
	d_v = ((1+z1)**2 * d_a**2 * z1 / h)**(1./3)

	pk_dirname='../../../demo_output_1/matter_power_lin/'
	k = np.loadtxt(pk_dirname+"k_h.txt")
	z2 = np.loadtxt(pk_dirname+"z.txt")
	p = np.loadtxt(pk_dirname+"p_k.txt")

	H0 = 70.0

	common = "data/wigglez_nov11_common.dataset"
	data_file = "data/wigglez_nov11a.dataset"
	root = '.'
	dataset = WigglezDataset(common, data_file, root)
	print dataset.likelihood(H0, z1[::-1], d_v[::-1], k, z2, p)


if __name__ == '__main__':
	test()

