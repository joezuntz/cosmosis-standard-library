from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
from astropy.io import fits
from scipy.interpolate import interp1d
import numpy as np
import os
import sys

# We need the modules below from the parent directory
dirname = os.path.split(__file__)[0]
sys.path.append(os.path.join(dirname, os.pardir))


from twopoint_cosmosis import theory_names, type_table
import twopoint
import gaussian_covariance


# Kids, don't name your python modules starting with a number.
# One day you will need to import them
twopt_like = __import__("2pt_like")


class TwoPointGammatMargLikelihood(twopt_like.TwoPointLikelihood):
    """
    This subclass of the 2pt likelihood adds marginalization over
    point-mass contributions to the small-scale galaxy galaxy-lensing
    """
    like_name='2pt'
    def __init__(self, options):
        super(TwoPointGammatMargLikelihood, self).__init__(options)
        
    def build_data(self):
        filename = self.options.get_string('data_file')

        covmat_name = self.options.get_string("covmat_name", "COVMAT")
        
        # Suffixes to added on to two point data from e.g. different experiments
        suffix_string = self.options.get_string('suffixes', default="")
        if suffix_string == "":
            # If there are no suffixes provided, then we create a list of a single empty suffix                               
            suffixes = [""]
        else:
            suffixes_temp = suffix_string.split()
            suffixes = []
            for suffix in suffixes_temp:
                if (suffix.lower() == 'none'):
                    suffixes.append('')
                else:
                    suffixes.append("_" + suffix)
        self.suffixes = suffixes
        
        #Prior width for analytically marginalized parameters
        #If <1. assume infinite prior
        self.sigma_a = self.options.get_double( "sigma_a", -1. )

        #do point-mass marginalization
        self.do_pm_marg = self.options.get_bool("do_pm_marg", False)

        #marginalize over small-scale gamma_t shear ratios
        self.do_smallscale_marg = self.options.get_bool("do_smallscale_marg", False)

        #name of gammat extension in datavector file
        self.gammat_name = self.options.get_string("gammat_name", "gammat")

        #use sigma_crit_inv factors in point-mass marginalization
        self.do_pm_sigcritinv = self.options.get_bool("do_pm_sigcritinv", False)

        #don't apply det factor in likelihood
        self.no_det_fac = self.options.get_bool("no_det_fac", False)

        print("Doing point-mass marginalization:", self.do_pm_marg)
        print("Using sigma_crit_inv factors in pm-marg:", self.do_pm_sigcritinv)
        print("Doing small-scale marginalization:", self.do_smallscale_marg)

        # This is the main work - read data in from the file
        self.two_point_data = twopoint.TwoPointFile.from_fits(
            filename, covmat_name)

        self.rescale_area = self.options.get_double( "rescale_area", -1.)
        if self.rescale_area > 0.:
            print("rescaling covariance by factor 1./rescale_area = %f"%(1./self.rescale_area))
            self.two_point_data.covmat /= self.rescale_area

        # All the names of two-points measurements that were found in the data
        # file
        all_names = [spectrum.name for spectrum in self.two_point_data.spectra]

        # We may not want to use all the likelihoods in the file.
        # We can set an option to only use some of them
        data_sets = self.options.get_string("data_sets", default="all")
        if data_sets != "all":
            data_sets = data_sets.split()
            self.two_point_data.choose_data_sets(data_sets)

        # The ones we actually used.
        self.used_names = [
            spectrum.name for spectrum in self.two_point_data.spectra]

        # Check for scale cuts. In general, this is a minimum and maximum angle for
        # each spectrum, for each redshift bin combination. Which is clearly a massive pain...
        # but what can you do?
        scale_cuts = {}
        gammat_scale_cuts_orig = {}
        for name in self.used_names:
            s = self.two_point_data.get_spectrum(name)
            for b1, b2 in s.bin_pairs:
                option_name = "angle_range_{}_{}_{}".format(name, b1, b2)
                if self.options.has_value(option_name):
                    r = self.options.get_double_array_1d(option_name)
                    if (name == self.gammat_name) and self.do_smallscale_marg:
                        #save gammat scale cuts separately in case doing 
                        gammat_scale_cuts_orig[(b1,b2)] = r
                    else:
                        scale_cuts[(name, b1, b2)] = r

        #if we doing the small scale shear ratio shizz
        #then the minimum scale cut needs to be the same
        #for all lens-source pairs for a given lens bin
        self.gammat_scale_cuts = {}
        if self.do_smallscale_marg:
            for (lens_bin, source_bin),scale_cut in gammat_scale_cuts_orig.items():
                if lens_bin not in self.gammat_scale_cuts:
                    self.gammat_scale_cuts[lens_bin] = scale_cut
                else:
                    if scale_cut != self.gammat_scale_cuts[lens_bin]:
                        print("ignoring gammat minimum scale cut for bin pair %d,%d = %f"%(lens_bin, source_bin, scale_cut))
                        print("using %f instead"%(self.gammat_scale_cuts[lens_bin]))

        # Now check for completely cut bins
        # example:
        # cut_wtheta = 1,2  1,3  2,3
        bin_cuts = []
        for name in self.used_names:
            s = self.two_point_data.get_spectrum(name)
            option_name = "cut_{}".format(name)
            if self.options.has_value(option_name):
                cuts = self.options[option_name].split()
                cuts = [eval(cut) for cut in cuts]
                for b1, b2 in cuts:
                    bin_cuts.append((name, b1, b2))

        if scale_cuts or bin_cuts:
            self.two_point_data.mask_scales(scale_cuts, bin_cuts)
        else:
            print("No scale cuts mentioned in ini file.")

        # Info on which likelihoods we do and do not use
        print("Found these data sets in the file:")
        total_data_points = 0
        final_names = [
            spectrum.name for spectrum in self.two_point_data.spectra]
        for name in all_names:
            if name in final_names:
                data_points = len(self.two_point_data.get_spectrum(name))
            else:
                data_points = 0
            if name in self.used_names:
                print("    - {}  {} data points after cuts {}".format(
                    name,  data_points, "  [using in likelihood]"))
                total_data_points += data_points
            else:
                print("    - {}  {} data points after cuts {}".format(
                    name, data_points, "  [not using in likelihood]"))
        print("Total data points used = {}".format(total_data_points))

        # Convert all units to radians.  The units in cosmosis are all
        # in radians, so this is the easiest way to compare them.
        for spectrum in self.two_point_data.spectra:
            if spectrum.is_real_space():
                spectrum.convert_angular_units("rad")

        # build up the data vector from all the separate vectors.
        # Just concatenation
        datavector = np.concatenate(
            [spectrum.value for spectrum in self.two_point_data.spectra])

        # Make sure
        if len(datavector) == 0:
            raise ValueError(
                """No data was chosen to be used from 2-point data file {0}. 
                It was either not selected in data_sets or cut out""".format(filename))

        #get covariance - this is what we start from before any analytic marginalization
        self.cov_orig = self.build_covariance()
        self.inv_cov_orig = np.linalg.inv( self.cov_orig )
        print("Original covariance log(det):", np.linalg.slogdet(self.cov_orig))
        self.logdet_fac = 0.
        #use self.det_fac to record changes in covariance determinant when we've applied
        #analytic marginalization terms - this is required to normalize the likelihood
        #number of points in the datavector
        self.n_datavector = len(datavector)

        #If we're doing analytic marginalization, we need to generate the template matrices
        #If we're not using sigma_crit_inv factors, we can generate the new inverse covariance here
        if self.do_pm_marg or self.do_smallscale_marg:
            if self.gammat_name in self.used_names:
                #find the datavector index, gammat_start_ind, at which gammat starts
                gammat_start_ind=0
                for s in self.two_point_data.spectra:
                    if s.name != self.gammat_name:
                        gammat_start_ind += len(s.value)
                    else:
                        gammat_spec = s
                        break

                #if we're retaining shear ratio info we need the sigma_crit_inv factors
                #we can't read them here since they're cosmology-dependent, but do get the section name
                #for later use
                if self.do_pm_sigcritinv or self.do_smallscale_marg:
                    lens_sample = (gammat_spec.kernel1).replace("nz_","")
                    source_sample = (gammat_spec.kernel2).replace("nz_","")
                    self.sigma_crit_inv_section = "sigma_crit_inv_%s_%s"%(lens_sample, source_sample)
                    print("will get sigma_crit_inv factors from section %s"%(self.sigma_crit_inv_section))
                    #also if we're using sigma_crit_inv factors then the covariance is cosmology-dependent.
                    self.constant_covariance = False
                    print("using sigma_crit_inv factors to retain shear-ratio info")

                #Now loop through lens-source pairs constructing the template matrices.
                bin_pairs = gammat_spec.get_bin_pairs()
                self.lens_bin_ids = list(set([p[0] for p in bin_pairs]))
                self.source_bin_ids = list(set([p[1] for p in bin_pairs]))

                if self.do_pm_marg:
                    self.bin_pair_dv_inds = {}
                    if self.do_pm_sigcritinv:
                        #the template matrix is n_lens_bin x N_data 
                        template_matrix_pm = np.zeros( (len(self.lens_bin_ids), self.n_datavector ) )
                        print("setting up template matrix for pm-marginalization with shape:", template_matrix_pm.shape)
                        #loop through lens bins generating the rows of the template matrix
                        for i,lens_bin in enumerate(self.lens_bin_ids):
                            template_vector = np.zeros(self.n_datavector)
                            lens_use = (gammat_spec.bin1==lens_bin)
                            #if we're doing the small-scale marg, we don't want those scales in
                            #the pm template matrix
                            if self.do_smallscale_marg:
                                lens_use *= gammat_spec.angle > np.radians(self.gammat_scale_cuts[(lens_bin)][0]/60.)
                            lens_spec_inds, = np.where(lens_use)
                            lens_dv_inds = lens_spec_inds + gammat_start_ind
                            theta_m2 = gammat_spec.angle[lens_use]**-2
                            template_vector[lens_dv_inds] = theta_m2
                            #Now loop through source bins getting bin_pair_dv_inds to get the sigma_crit_inv factors later
                            for source_bin in self.source_bin_ids:
                                bin_pair_use = lens_use * (s.bin2==source_bin)
                                bin_pair_spec_inds, = np.where(bin_pair_use)
                                self.bin_pair_dv_inds[lens_bin, source_bin] = bin_pair_spec_inds + gammat_start_ind
                            template_matrix_pm[i] = template_vector
                    
                        self.template_matrix_pm = template_matrix_pm

                    else:
                        template_matrix = np.zeros( ( len(bin_pairs), self.n_datavector ) )
                        print("setting up template matrix for pm-marg with shape:", template_matrix.shape)
                        for ibin_pair, bin_pair in enumerate(bin_pairs):
                            lens_bin,source_bin = bin_pair
                            use = (gammat_spec.bin1==lens_bin)*(gammat_spec.bin2==source_bin)
                            if self.do_smallscale_marg:
                                use *= gammat_spec.angle > np.radians(self.gammat_scale_cuts[(lens_bin)][0]/60.)
                            spec_inds, = np.where(use)
                            dv_inds = spec_inds + gammat_start_ind
                            theta_vals = gammat_spec.angle[use]
                            theta_m2 = theta_vals**-2
                            template_vector = theta_m2
                            template_matrix[ibin_pair][dv_inds] = template_vector
                        self.template_matrix_pm = template_matrix

                        if not self.do_smallscale_marg:
                            # We're not using sigma_crit_inv factors for the point-mass 
                            # marginalization, and we're not doing the small-scale marginalization
                            # so nothing cosmology-dependent, so we can go ahead and doctor the covariance now
                            # using  Woodbury matrix identity
                            U, V = template_matrix.T, template_matrix
                            UCinvV = np.matmul(V, np.matmul(self.inv_cov_orig, U))
                            if self.sigma_a>0.:
                                print("using sigma_a = %f"%self.sigma_a)
                                X = self.sigma_a**-2 * np.identity(UCinvV.shape[0]) + UCinvV 
                            else:
                                print("taking limit sigma_a -> inf")
                                X = UCinvV
                            Xinv = np.linalg.inv(X)

                            Y = np.matmul( Xinv, np.matmul( V, self.inv_cov_orig ) )
                            sub_from_inv = np.matmul( self.inv_cov_orig, np.matmul( U, Y ))
                            print('subtracting from inverse covariance:')
                            print(sub_from_inv)
                            self.inv_cov = self.inv_cov_orig - sub_from_inv

                            #det factor 
                            s, logdet_fac = np.linalg.slogdet(X)
                            assert s>0
                            self.logdet_fac = logdet_fac
                            if self.sigma_a > 0:
                                self.logdet_fac += np.log(self.sigma_a**2)
                            print("logdet_fac: %f"%self.logdet_fac)

                else:
                    self.template_matrix_pm = None

                #Now setup small-scale marginalization
                if self.do_smallscale_marg:
                    self.template_matrix_smallscale_info = []
                    #loop through lens bins and angular bins for that lens bin
                    #collect datavector indices and redshift bin pair for each angular bin being cut
                    #this is the info required to build the template matrix during execute
                    unq_angular_bins, unq_angular_bin_inds = np.unique( gammat_spec.angular_bin, return_index=True )
                    unq_angle_vals = gammat_spec.angle[ unq_angular_bin_inds ]

                    for lens_bin in self.lens_bin_ids:
                        for angular_bin, angle_val in zip(unq_angular_bins, unq_angle_vals):
                            if angle_val > np.radians(self.gammat_scale_cuts[(lens_bin)][0]/60.):
                                continue
                            else:
                                template_vector = np.zeros(self.n_datavector)
                            gammat_spec_use = (gammat_spec.bin1==lens_bin)*(gammat_spec.angular_bin==angular_bin)
                            gammat_spec_inds, = np.where(gammat_spec_use)
                            source_bins = gammat_spec.bin2[gammat_spec_use]
                            dv_inds = gammat_start_ind + gammat_spec_inds
                            self.template_matrix_smallscale_info.append((dv_inds, lens_bin, source_bins))

                    print("set up template matrix info for small-scale marginalization with shape: (%d,%d)"%(len(self.template_matrix_smallscale_info), self.n_datavector) )
                else:
                    self.template_matrix_smallscale_info = None
            else:
                print("No %s in datavector so not doing any analytic marginalization"%self.gammat_name)
                self.inv_cov = np.linalg.inv(self.cov_orig)
        else:
            print("Not doing any analytic marginalization")
            self.inv_cov = np.linalg.inv(self.cov_orig)

        return None, datavector


    def extract_covariance(self, block):
        self.chol = np.linalg.cholesky(self.cov_orig)
        return self.cov_orig

    def extract_inverse_covariance(self, block):
        if (self.do_smallscale_marg) or (self.do_pm_sigcritinv and self.do_pm_marg):
            #We've constructed the template matrix in build_data as part of the module setup
            #Now we need to loop through source bins applying sigma crit inv factors - have to
            #do this now since these are cosmology dependent. 
            #read in sigma_crit_inv factors
            #set self.log_det_fac to zero to start with, it can pick up factors from both the pm and small-scale
            #marginalization.
            self.logdet_fac = 0.
            sigma_crit_inv_dict = {}
            for lens_bin in self.lens_bin_ids:
                for source_bin in self.source_bin_ids:
                    sigma_crit_inv_dict[lens_bin,source_bin] = block.get_double(self.sigma_crit_inv_section,
                                                            "sigma_crit_inv_%d_%d"%(lens_bin, source_bin))
            #apply sigma crit inv factors for pm-marg template matrix
            if self.do_pm_marg:
                template_matrix_pm = np.copy(self.template_matrix_pm)
                if self.do_pm_sigcritinv:
                    for i, lens_bin in enumerate(self.lens_bin_ids):
                        for j, source_bin in enumerate(self.source_bin_ids):
                            sig_crit_inv = sigma_crit_inv_dict[lens_bin,source_bin]
                            #template_inds = (self.source_bin_vector == source_bin) * self.lens_use_list[i]
                            inds = self.bin_pair_dv_inds[lens_bin, source_bin]
                            template_matrix_pm[i][inds] *= sig_crit_inv 
            else:
                template_matrix_pm = None

            if self.do_smallscale_marg:
                template_matrix_smallscale = np.zeros( (len(self.template_matrix_smallscale_info), self.n_datavector) )
                for i,(dv_inds, lens_bin, source_bins) in  enumerate(self.template_matrix_smallscale_info):
                    for dv_ind, source_bin in zip(dv_inds, source_bins):
                        template_matrix_smallscale[i][dv_ind] = sigma_crit_inv_dict[lens_bin, source_bin]
                if template_matrix_pm is not None:
                    template_matrix_with_sig = np.vstack((template_matrix_pm, template_matrix_smallscale))
                else:
                    template_matrix_with_sig = template_matrix_smallscale
            else:
                template_matrix_with_sig = template_matrix_pm
            
            #Now construct inverse covariance
            U, V = template_matrix_with_sig.T, template_matrix_with_sig
            UCinvV = np.matmul(V, np.matmul(self.inv_cov_orig, U))
            if self.sigma_a>0.:
                X = self.sigma_a**-2 * np.identity(UCinvV.shape[0]) + UCinvV 
            else:
                #infinite prior case
                X = UCinvV
            Xinv = np.linalg.inv(X)
            Y = np.matmul( Xinv, np.matmul( V, self.inv_cov_orig ) )
            sub_from_inv = np.matmul( self.inv_cov_orig, np.matmul( U, Y ))
            self.inv_cov = self.inv_cov_orig - sub_from_inv
            #we're now changing the covariance in a cosmology dependent way...
            #The determinant is changed by a factor det(I + U^TC^{-1}V )
            s, logdet_fac = np.linalg.slogdet(X)
            try:
                assert s>0
            except AssertionError as e:
                print(s, logdet_fac)
                print(X)
                raise(e)
            self.logdet_fac += logdet_fac
            if self.sigma_a > 0:
                self.logdet_fac += np.log(self.sigma_a ** 2)

        if self.no_det_fac:
            self.logdet_fac = 0.

        return self.inv_cov

    def extract_covariance_log_determinant(self, block):
        sign, log_det = np.linalg.slogdet(self.cov_orig)
        log_det += self.logdet_fac
        return log_det

    def build_inverse_covariance(self):
        """
        Override the build_inverse_covariance method to change
        how the inverse is generated from the covariance.

        When the covariance is generated from a suite of simulations,
        for example, the simple inverse is not the best estimate.

        """
        return self.inv_cov

setup, execute, cleanup = TwoPointGammatMargLikelihood.build_module()
