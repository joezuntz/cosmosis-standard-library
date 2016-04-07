"""
This is a general 2-point likelihood designed to connect
cosmosis to the fits data format developed for the Dark
Energy Survey 2-point data.

At the moment it is designed to handle C_ell likelihoods 
for shear, galaxy-galaxy lensing, and clustering in 
photometric surveys, but its behaviour is general enough
(with the possible exception of the covariance code) to
be extended to other two-point analyses.

"""

from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
from scipy.interpolate import interp1d
import numpy as np
import twopoint
import gaussian_covariance
from twopoint_cosmosis import theory_names, type_table

default_array = np.repeat(-1.0, 99)
def is_default(x):
	return len(x)==len(default_array) and (x==default_array).all()

def convert_nz_steradian(n):
    return n * (41253.0*60.*60.) / (4*np.pi)



class TwoPointLikelihood(GaussianLikelihood):
	like_name = "2pt"

	def __init__(self, options):
		#We may decide to use an analytic gaussian covariance
		#in that case we won't load the covmat.
		self.gaussian_covariance = options.get_bool("gaussian_covariance", False)
		if self.gaussian_covariance:
			self.constant_covariance = False

		super(TwoPointLikelihood, self).__init__(options)

	
	def build_data(self):
		filename = self.options.get_string('data_file')
		self.do_plot = self.options.get_bool('do_plot', default=False)

		if self.gaussian_covariance:
			covmat_name = None
			area = self.options.get_double("survey_area") # in square degrees
			self.sky_area = area * (np.pi*np.pi)/(180*180)

			def get_arr(x):
				if self.options.has_value(x):
					a = self.options[x]
					if not isinstance(a, np.ndarray):
						a = [a]
				else:
					a = default_array
				return a

			self.number_density_shear_bin = get_arr("number_density_shear_bin")
			self.number_density_lss_bin = get_arr("number_density_lss_bin")
			self.sigma_e_bin = get_arr("sigma_e_bin")

		else:
			covmat_name = self.options.get_string("covmat_name", "COVMAT")

		# This is the main work - read data in from the file
		self.two_point_data = twopoint.TwoPointFile.from_fits(filename, covmat_name)

		#Potentially cut out lines. For some reason one version of 
		#this file used zeros to mark masked values.
		if self.options.get_bool("cut_zeros", default=False):
			print "Removing 2-point values with value=0.0"
			self.two_point_data.mask_bad(0.0)


		#All the names of two-points measurements that were found in the data
		#file
		all_names = [spectrum.name for spectrum in self.two_point_data.spectra]

		#We may not want to use all the likelihoods in the file. 
		#We can set an option to only use some of them
		data_sets = self.options.get_string("data_sets", default="all")
		if data_sets!="all":
			data_sets = data_sets.split()
			self.two_point_data.choose_data_sets(data_sets)

		#The ones we actually used. 
		used_names = [spectrum.name for spectrum in self.two_point_data.spectra]

		ell_max = self.options.get_double("ell_max", default=-1.0)
		spectra_to_cut = self.options.get_string("spectra_to_cut", default="all")
		if spectra_to_cut != "all":
			spectra_to_cut = spectra_to_cut.split()
		if ell_max>=0:
			self.two_point_data.mask_scale(spectra_to_cut, max_scale=ell_max)

		#Info on which likelihoods we do and do not use
		print "Found these data sets in the file:"
		for name in all_names:
			if name in used_names:
				print "    - ",name, "  [using in likelihood]"
			else:
				print "    - ",name, "  [not using in likelihood]"

		#build up the data vector from all the separate vectors.
		#Just concatenation
		data_vector = np.concatenate([spectrum.value for spectrum in self.two_point_data.spectra])

		#Make sure
		if len(data_vector)==0:
			raise ValueError("No data was chosen to be used from 2-point data file {0}. It was either not selectedin data_sets or cut out".format(filename))

		#The x data is not especially useful here, so return None.
		#We will access the self.two_point_data directly later to 
		#determine ell/theta values
		return None, data_vector


	def build_covariance(self):
		C = np.array(self.two_point_data.covmat)		
		r = self.options.get_int('covariance_realizations', default=0)
		if r:
			p = C.shape[0]
			print "You set covariance_realizations={} in the 2pt likelihood parameter file".format(r)
			print "So I will apply the Anderson-Hartlap correction to the covariance matrix"
			print "The covariance matrix is {}x{}".format(p,p)
			#This x is the inverse of the alpha used in the old code
			#because that applied to the weight matrix not the covariance
			x = (r - 1.0) / (r - p - 2.0)
			print "So the correction scales the covariance matrix by (r - 1) / (r - n - 2) = {}".format(x)
			C = C * x

		return C
		

	def extract_theory_points(self, block):
		theory = []
		#We may want to save these splines for the covariance matrix later
		self.theory_splines = {}


		# We have a collection of data vectors, one for each spectrum
		# that we include. We concatenate them all into one long vector,
		# so we do the same for our theory data so that they match

		#We will also save angles and bin indices for plotting convenience,
		#although these are not actually used in the likelihood
		angle = []
		bin1 = []
		bin2 = []

		#Now we actually loop through our data sets
		for spectrum in self.two_point_data.spectra:
			theory_vector, angle_vector, bin1_vector, bin2_vector = self.extract_spectrum_prediction(block, spectrum)
			theory.append(theory_vector)
			angle.append(angle_vector)
			bin1.append(bin1_vector)
			bin2.append(bin2_vector)

		#We also collect the ell or theta values.
		#The gaussian likelihood code itself is not expecting these,
		#so we just save them here for convenience.
		angle = np.concatenate(angle)
		bin1 = np.concatenate(bin1)
		bin2 = np.concatenate(bin2)
		block[names.data_vector, self.like_name+"_angle"] = angle
		block[names.data_vector, self.like_name+"_bin1"] = bin1
		block[names.data_vector, self.like_name+"_bin2"] = bin2

		#the thing it does want is the theory vector, for comparison with
		#the data vector
		theory = np.concatenate(theory)
		return theory

	def extract_spectrum_prediction(self, block, spectrum):
		#We may need theory predictions for multiple different
		#types of spectra: e.g. shear-shear, pos-pos, shear-pos.
		#So first we find out from the spectrum where in the data
		#block we expect to find these - mapping spectrum types
		#to block names
		section, x_name, y_name = theory_names(spectrum)

		#We need the angle (ell or theta depending on the spectrum)
		#for the theory spline points - we will be interpolating
		#between these to get the data points
		angle_theory = block[section, x_name]

		#Now loop through the data points that we have.
		#For each one we have a pairs of bins and an angular value.
		#This assumes that we can take a single sample point from
		#each theory vector rather than integrating with a window function
		#over the theory to get the data prediction - this will need updating soon.
		bin_data = {}
		theory_vector = []

		#For convenience we will also return the bin and angle (ell or theta)
		#vectors for this bin too.
		angle_vector = []
		bin1_vector = []
		bin2_vector = []
		for (b1, b2, angle) in zip(spectrum.bin1,spectrum.bin2,spectrum.angle):
			#We are going to be making splines for each pair of values that we need.
			#We make splines of these and cache them so we don't re-make them for every
			#different theta/ell data point
			if (b1,b2) in bin_data:
				#either use the cached spline
				theory_spline = bin_data[(b1,b2)]
			else:
				#or make a new cache value
				#load from the data block and make a spline
				#and save it
				if block.has_value(section, y_name.format(b1,b2)):
					theory = block[section, y_name.format(b1,b2)]
				#It is okay to swap if the spectrum types are the same - symmetrical
				elif block.has_value(section, y_name.format(b2,b1)) and spectrum.type1==spectrum.type2:
					theory = block[section, y_name.format(b2,b1)]
				else:
					raise ValueError("Could not find theory prediction {} in section {}".format(y_name.format(b1,b2),section))
				theory_spline = interp1d(angle_theory, theory)
				bin_data[(b1,b2)] = theory_spline
				#This is a bit silly, and is a hack because the
				#book-keeping is very hard.
				bin_data[y_name.format(b1,b2)] = theory_spline

			#use our spline - interpolate to this ell or theta value
			#and add to our list
			theory = theory_spline(angle)
			theory_vector.append(theory)
			angle_vector.append(angle)
			bin1_vector.append(b1)
			bin2_vector.append(b2)

		#We are saving the theory splines as we may need them
		#to calculate covariances later
		self.theory_splines[section] = bin_data

		if self.do_plot:
			import pylab
			nbin = max(spectrum.nbin(), spectrum.nbin())
			for b1 in xrange(1,nbin+1):
				for b2 in xrange(1,nbin+1):
					if (b1,b2) not in bin_data:
						continue
					pylab.subplot(nbin, nbin, (b1-1)*nbin+b2)
					pylab.loglog(angle_theory, bin_data[(b1,b2)](angle_theory))
					xdata, ydata = spectrum.get_pair(b1, b2)
					yerr = spectrum.get_error(b1, b2)
					pylab.errorbar(xdata, ydata, yerr, fmt='+')
					pylab.xlim(xdata.min(), xdata.max())
					pylab.ylim(ydata.min(), ydata.max())
			pylab.savefig("cmp_{}.pdf".format(spectrum.name))
			pylab.close()

			
		#Return the whole collection as an array
		theory_vector = np.array(theory_vector)

		#For convenience we also save the angle vector (ell or theta)
		#and bin indices
		angle_vector = np.array(angle_vector)
		bin1_vector = np.array(bin1_vector, dtype=int)
		bin2_vector = np.array(bin2_vector, dtype=int)

		return theory_vector, angle_vector, bin1_vector, bin2_vector

	def extract_covariance(self, block):
		assert self.gaussian_covariance, "Set constant_covariance=F but somehow not with Gaussian covariance.  Internal error - please open an issue on the cosmosis site."
		

		C = []
		# s and t index the spectra that we have. e.g. s or t=1 might be the full set of 
		#shear-shear measuremnts
		for s,AB in enumerate(self.two_point_data.spectra[:]):
			M = []
			for t,CD in enumerate(self.two_point_data.spectra[:]):
				print "Looking at covariance between {} and {} (s={}, t={})".format(AB.name, CD.name, s, t)
				#We only calculate the upper triangular.
				#Get the lower triangular here. We have to 
				#transpose it compared to the upper one.
				if s>t:
					MI = C[t][s].T
				else:
					MI = gaussian_covariance.compute_gaussian_covariance(self.sky_area, 
						self._lookup_theory_cl, block, AB, CD)
				M.append(MI)
			C.append(M)

		#C is now a list of lists of 2D arrays.
		#Now turn C into a big 2D array by stacking
		#the arrays
		C = np.vstack([np.hstack(CI) for CI in C])

		return C

	def _lookup_theory_cl(self, block, A, B, i, j, ell):
		"""
		This is a helper function for the compute_gaussian_covariance code.
		It looks up the theory value of C^{ij}_{AB}(ell) in the 
		"""
		# We have already saved splines into the theory space earlier
		# when constructing the theory vector.
		# So now we just need to look those up again, using the same
		# code we use in the twopoint library.
		section, ell_name, value_name = type_table[A, B]
		assert ell_name=="ell", "Gaussian covariances are currently only written for C_ell, not other 2pt functions"
		d = self.theory_splines[section]

		#We save the splines with these names when we extract the theory vector
		name_ij = value_name.format(i,j)
		name_ji = value_name.format(j,i)

		#Hopefully we already have the theory spline extracted
		if name_ij in d:
			spline = d[name_ij]
		#For symmetric spectra (not just auto-correlations, but any thing like C_EE or C_NN where
		#we cross-correlate something with itself) we can use ji for ij as it is the same. This is 
		#not true for cross spectra
		elif name_ji in d and (A==B):
			spline = d[name_ji]
		else:
		#It's possible too that we need something for the covariance that we didn't need for the 
		#data vector - for example to got the covariance between C^EE and C^NN we need C^NE even
		#if we don't have any actual measurements of NE. In that case we have to g
			angle_theory = block[section, ell_name]
			if block.has_value(section, name_ij):
				theory = block[section, name_ij]
			# The same symmetry argument as above applies
			elif block.has_value(section, name_ji) and A==B:
				theory = block[section, name_ji]
			else:
				raise ValueError("Could not find theory prediction {} in section {}".format(value_name.format(i,j), section))

			spline = interp1d(angle_theory, theory)
			#Finally cache this so we don't have to do this again.
			d[name_ij] = spline

		obs_cl = spline(ell)

		#For shear-shear the noise component is sigma^2 / number_density_bin
		#and for position-position it is just 1/number_density_bin
		if (A==B) and (A==twopoint.Types.galaxy_shear_emode_fourier.name) and (i==j):
			if i>len(self.number_density_shear_bin) or i>len(self.sigma_e_bin) or is_default(self.sigma_e_bin) or is_default(self.number_density_shear_bin):
				raise ValueError("Not enough number density bins for shear specified")
			noise = self.sigma_e_bin[i-1]**2 / convert_nz_steradian(self.number_density_shear_bin[i-1])
			obs_cl += noise
		if (A==B) and (A==twopoint.Types.galaxy_position_fourier.name) and (i==j):
			if i>len(self.number_density_lss_bin) or is_default(self.number_density_lss_bin):
				raise ValueError("Not enough number density bins for lss specified")
			noise = 1.0 / convert_nz_steradian(self.number_density_lss_bin[i-1])
			obs_cl += noise


		return obs_cl
		




setup, execute, cleanup = TwoPointLikelihood.build_module()
