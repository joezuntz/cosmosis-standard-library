#coding:utf-8
import os
import ctypes as ct
import numpy as np
import limber
from gsl_wrappers import GSLSpline, NullSplineError, GSLSpline2d, BICUBIC
from cosmosis.datablock import names, option_section, BlockError
from enum34 import Enum
import re
import sys
import scipy.interpolate as interp

class Power3D(object):
    """
    Class representing the 3D power spectrum that enters the Limber calculation.
    Most Power spectra are source-specific, like intrinsic alignments and galaxy
    density.  Others are generic, like matter_power (linear and nonlinear)

    """
    # Must specify these:
    # section = "?????"
    # source_specific = None
    def __init__(self, suffix=""):
        if self.source_specific:
            self.section_name = self.section + suffix
        else:
            self.section_name = self.section
    def __hash__(self):
        return hash(self.section_name)

class MatterPower3D(Power3D):
    section = "matter_power_nl"
    source_specific = False

class LinearMatterPower3D(Power3D):
    section = "matter_power"
    source_specific = False

class GalaxyPower3D(Power3D):
    section = "galaxy_power"
    source_specific = True

class IntrinsicPower3D(Power3D):
    section = "intrinsic_power"
    source_specific = True

class IntrinsicBBPower3D(Power3D):
    section = "intrinsic_power_bb"
    source_specific = True

class MatterGalaxyPower3D(Power3D):
    section = "matter_galaxy_power"
    source_specific = True

class MatterIntrinsicPower3D(Power3D):
    section = "matter_intrinsic_power"
    source_specific = True

class GalaxyIntrinsicPower3D(Power3D):
    section = "galaxy_intrinsic_power"
    source_specific = True


class Spectrum(object):
    autocorrelation = False
    #These should make it more obvious if the values are not overwritten by subclasses
    power_3d_type = "?"
    kernels = "??"
    name = "?"
    prefactor_power = np.nan
    #the default is no magnification. If the two fields are magnification terms
    #that should pick up factors of 2 alpha_i - 1
    #then subclasses should include 1 and/or 2 in this.
    magnification_prefactors = {}

    def __init__(self, source, power_3d, sample_a=names.wl_number_density, 
                 sample_b=names.wl_number_density, save_name=""):
        #caches of n(z), w(z), P(k,z), etc.
        self.source = source
        self.power_3d = power_3d
        self.sample_a = sample_a
        self.sample_b = sample_b
        self.save_name=save_name

    def get_name(self):
        if self.save_name:
            return self.name+"_"+self.save_name
        return self.name

    def nbins(self):
        na = len(self.source.kernels_A[self.kernels[ 0]+"_"+self.sample_a])
        nb = len(self.source.kernels_B[self.kernels[-1]+"_"+self.sample_b])
        return na,nb

    def is_autocorrelation(self):
        """
        This is an autocorrelation if the basic type is an auto-correlation
        (e.g. shear-shear, position-position, but not shear position)
        and the two n(z) samples are the same.
        """
        return self.autocorrelation and (self.sample_a==self.sample_b)

    @classmethod
    def option_name(cls):
        """Convert the CamelCase name to hypen-separated.
        For example ShearShear becomes shear-shear
        """
        name = cls.__name__
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1-\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1-\2', s1).lower()

    def prefactor(self, block, bin1, bin2):
        if self.prefactor_power == 0:
            #This should be okay for the magnification terms
            #since none of those have prefactor_power==0
            return 1.0
        c_kms = 299792.4580
        omega_m = block[names.cosmological_parameters, "omega_m"]
        shear_scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m
        f = shear_scaling**self.prefactor_power
        if 1 in self.magnification_prefactors:
            f *= self.magnification_prefactor(block, bin1)
        if 2 in self.magnification_prefactors:
            f *= self.magnification_prefactor(block, bin2)
        return f

    #Some spectra include a magnification term for one or both bins.
    #In those cases an additional term from the galaxy luminosity
    #function is included in the prefactor. 
    def magnification_prefactor(self, block, bin):
        #Need to move this so that it references a particular
        #galaxy sample not the generic galaxy_luminosity_function
        alpha = np.atleast_1d(np.array(block[names.galaxy_luminosity_function, "alpha_binned"]))
        return 2*(alpha[bin]-1)


    def compute(self, block, ell, bin1, bin2, relative_tolerance=1.e-3, 
                absolute_tolerance=0., do_extended_limber=False, term=None):
        #Get the required kernels
        #maybe get from rescaled version of existing spectra
        P, D = self.get_power_growth(block, bin1, bin2)
        K1 = self.source.kernels_A[self.kernels[ 0]+"_"+self.sample_a][bin1]
        K2 = self.source.kernels_B[self.kernels[-1]+"_"+self.sample_b][bin2]


        K = block.get_double(names.cosmological_parameters, "K")

        #Call the integral
        #if do_extended_limber is True, pass ext_order=1
        ext_order=0
        if do_extended_limber:
            ext_order=1
        c_ell = limber.extended_limber(K1, K2, P, D, ell, self.prefactor(block, bin1, bin2), 
            rel_tol=relative_tolerance, abs_tol=absolute_tolerance, ext_order=ext_order, K=K)
        self.clean_power_growth(P, D)
        return c_ell

    def get_power_growth(self, block, bin1, bin2):
        return self.source.power[self.power_3d], self.source.growth[self.power_3d]

    def clean_power_growth(self, P, D):
        #This gets done later for the base class
        return 0

    def kernel_peak(self, block, bin1, bin2, a_of_chi):
        K1 = self.source.kernels_A[self.kernels[ 0]+"_"+self.sample_a][bin1]
        K2 = self.source.kernels_B[self.kernels[-1]+"_"+self.sample_b][bin2]
        chi_peak = limber.get_kernel_peak(K1, K2)
        a_peak = a_of_chi(chi_peak)
        z_peak = 1./a_peak - 1
        return chi_peak, z_peak

    def prep_spectrum(self, *args, **kwargs):
        #no prep required for base class
        return 0

# This is pretty cool.
# You can make an enumeration class which
# contains a list of possible options for something.
# But the options can each be anything - full python objects, and in this
# case they are classes of spectrum. So we can easily look these up by name,
# loop through them, etc.
class SpectrumType(Enum):

    class ShearShear(Spectrum):
        power_3d_type = MatterPower3D
        kernels = "W W"
        autocorrelation = True
        name = names.shear_cl
        prefactor_power = 2

    class ShearIntrinsic(Spectrum):
        power_3d_type = MatterIntrinsicPower3D
        kernels = "W N"
        autocorrelation = False
        name = names.shear_cl_gi
        prefactor_power = 1

    class IntrinsicIntrinsic(Spectrum):
        power_3d_type = IntrinsicPower3D
        kernels = "N N"
        autocorrelation = True
        name = names.shear_cl_ii
        prefactor_power = 0

    class IntrinsicbIntrinsicb(Spectrum):
        power_3d_type = IntrinsicBBPower3D
        kernels = "N N"
        autocorrelation = True
        name = "shear_cl_bb"
        prefactor_power = 0

    class PositionPosition(Spectrum):
        power_3d_type = GalaxyPower3D
        kernels = "N N"
        autocorrelation = True
        name = "galaxy_cl"
        prefactor_power = 0

    class MagnificationPosition(Spectrum):
        power_3d_type = MatterGalaxyPower3D
        kernels = "W N"
        autocorrelation = False
        name = "magnification_galaxy_cl"
        prefactor_power = 1    
        magnification_prefactors = (1,)

    #
    class MagnificationMagnification(Spectrum):
        power_3d_type = MatterPower3D
        kernels = "W W"
        autocorrelation = True
        name = "magnification_cl"
        prefactor_power = 2
        magnification_prefactors = (1,2)


    class PositionShear(Spectrum):
        power_3d_type = MatterGalaxyPower3D
        kernels = "N W"
        autocorrelation = False
        name = "galaxy_shear_cl"
        prefactor_power = 1

    class PositionIntrinsic(Spectrum):
        power_3d_type = GalaxyIntrinsicPower3D
        kernels = "N N"
        autocorrelation = False
        name = "galaxy_intrinsic_cl"
        prefactor_power = 0

    class MagnificationIntrinsic(Spectrum):
        power_3d_type = MatterIntrinsicPower3D
        kernels = "W N"
        autocorrelation = False
        name = "magnification_intrinsic_cl"
        prefactor_power = 1
        magnification_prefactors = (1,)

    class MagnificationShear(Spectrum):
        power_3d_type = MatterPower3D
        kernels = "W W"
        autocorrelation = False
        name = "magnification_shear_cl"
        prefactor_power = 2
        magnification_prefactors = (1,)

    class ShearCmbkappa(Spectrum):
        power_3d_type = MatterPower3D
        kernels = "W K"
        autocorrelation = False
        name = "shear_cmbkappa_cl"
        prefactor_power = 2

    class CmbkappaCmbkappa(Spectrum):
        power_3d_type = MatterPower3D
        kernels = "K K"
        autocorrelation = True
        name = "cmbkappa_cl"
        prefactor_power = 2

    class IntrinsicCmbkappa(Spectrum):
        power_3d_type = MatterIntrinsicPower3D
        kernels = "N K"
        autocorrelation = False
        name = "intrinsic_cmbkappa_cl"
        prefactor_power = 1

    class PositionCmbkappa(Spectrum):
        power_3d_type = MatterGalaxyPower3D
        kernels = "N K"
        autocorrelation = False
        name = "galaxy_cmbkappa_cl"
        prefactor_power = 1

class SpectrumCalculator(object):
    # It is useful to put this here so we can subclass to add new spectrum
    # types, for example ones done with modified gravity changes.
    spectrumType = SpectrumType
    def __init__(self, options):
        #General options
        self.verbose = options.get_bool(option_section, "verbose", False)
        self.fatal_errors = options.get_bool(option_section, "fatal_errors", False)
        self.get_kernel_peaks = options.get_bool(option_section, "get_kernel_peaks", False)
        self.save_kernel_zmax = options.get_double(option_section, "save_kernel_zmax", -1.0)
        self.do_extended_limber = options.get_bool(option_section, "do_extended_limber", False)
        self.do_rsd = options.get_bool(option_section, "do_rsd", False)
        self.nl_bias_section = options.get_string(option_section, "nl_bias_section", "bias_b1b2bs_perbin")
        #by default the nonlinearly biased spectra will start from the nonlinear matter power spectrum,
        #but in some situations you may want them to start from galaxy_power/galaxy_matter_power instead:
        self.use_galaxy_power = options.get_bool(option_section, "use_galaxy_power", False)

        # Check which spectra we are requested to calculate
        self.parse_requested_spectra(options)
        print("Will project these spectra into 2D:")
        for spectrum in self.req_spectra:
            print("    - ", spectrum.get_name())

        #Decide which kernels we will need to save.
        #The overall split is into A and B, the two samples we are correlating into.
        #Within each list of those we may need W(z) and N(z). These may be just
        #a single group for a single galaxy sample, or they may be more complicated.
        #The keys in kernels_A are things like "N_REDMAGIC" or "W_MAINSAMPLE", or if
        #everything is left with the defaults just "N_WL_NUMBER_DENSITY" and "W_WL_NUMBER_DENSITY"
        #The values are dictionaries where the keys are the bin numbers 1..nbin
        self.kernels_A = {}
        self.kernels_B = {}
        for spectrum in self.req_spectra:
            #names e.g. N_REDMAGIC, W_EUCLID
            kernel_a = spectrum.kernels[0] + "_" + spectrum.sample_a 
            kernel_b = spectrum.kernels[-1] + "_" + spectrum.sample_b 
            #one dictionary per named group
            self.kernels_A[kernel_a] = {}
            self.kernels_B[kernel_b] = {}


        self.power = {}
        self.growth = {}
        self.outputs = {}

        #And the req ell ranges.
        #We use log-spaced output
        ell_min = options.get_double(option_section, "ell_min")
        ell_max = options.get_double(option_section, "ell_max")
        n_ell = options.get_int(option_section, "n_ell")
        self.ell = np.logspace(np.log10(ell_min), np.log10(ell_max), n_ell)
        self.absolute_tolerance = options.get_double(option_section, "limber_abs_tol", 0.)
        self.relative_tolerance = options.get_double(option_section, "limber_rel_tol", 1.e-3)

    def parse_requested_spectra(self, options):
        #Get the list of spectra that we want to compute.
        #The full list
        self.req_spectra = []
        self.req_power = set()
 
        any_spectra_option_found = False
        for spectrum in self.spectrumType:
            spectrum = spectrum.value
            #By default we just do the shear-shear spectrum.
            #everything else is not done by default
            name = spectrum.option_name()
            try:
                value = options[option_section, name]
            #if value is not set at all, skip
            except BlockError:
                continue

            # If we find any of the options set then record that
            any_spectra_option_found = True

            # There are various ways a user can describe the spectra:
            # True/False, to use the default section names
            # string of form  euclid-lsst
            #  (one pair of n(z) samples to correlate,
            #   output in the default section)
            #  euclid-lsst:cross  euclid-euclid:auto
            #  (one or more pairs, output in the named section)

            if isinstance(value, bool):
                if value:
                    power3D = spectrum.power_3d_type()
                    self.req_spectra.append(spectrum(self, power3D))
                    print("Adding {}".format(power3D))
                    self.req_power.add(power3D)
                continue

            #Otherwise it must be a string - enforce this.
            if not (isinstance(value, str) or isinstance(value, unicode)):
                raise ValueError("Unknown form of value for option {} in project_2d: {}".format(name, value))

            value = value.strip()
            if not value:
                raise ValueError("Empty value for option {} in project_2d.".format(name))


            #now we are looking for things of the form
            #shear-shear = euclid-ska[:name]  
            #where we would now search for nz_euclid and nz_ska
            values = value.split()
            for value in values:
                try:
                    kernel_a, kernel_b = value.split('-',1)
                    #Optionally we can also name the spectrum, for example
                    # shear-shear = ska-ska:radio
                    # in which case the result will be saved into shear_cl_radio
                    # instead of just shear_cl.
                    # This will be necessary in the case that we run multiple spectra,
                    # e.g. 
                    # #shear-shear = euclid-ska:cross  euclid-euclid:optical  ska-ska:radio
                    # would be needed to avoid clashes. 
                    # Can also allow
                    # #intrinsic-intrinsic = des_source-des_source:des_power:des_cl kids_source-kids_source:kids_power:kids_cl
                    # to use the suffix XXX or YYY on the IA and galaxy density 3D power spectrum inputs
                    if ":" in kernel_b:
                        kernel_b, save_name=kernel_b.split(":",1)
                        if ":" in save_name:
                            power_suffix,save_name = save_name.split(':',1)
                            power_suffix = "_"+power_suffix
                        else:
                            power_suffix = ""
                    else:
                        save_name = ""
                        power_suffix = ""
                    kernel_a = kernel_a.strip()
                    kernel_b = kernel_b.strip()
                    #The self in the line below is not a mistake - the source objects
                    #for the spectrum class is the SpectrumCalculator itself
                    power3D = spectrum.power_3d_type(power_suffix)
                    self.req_power.add(power3D)
                    self.req_spectra.append(spectrum(self, power3D, kernel_a, kernel_b, save_name))
                    print("Calculating Limber: Kernel 1 = {}, Kernel 2 = {}, P_3D = {} --> Output: {}".format(
                        kernel_a, kernel_b, power3D.section_name, save_name))
                except:
                    raise
                    raise ValueError("To specify a P(k)->C_ell projection with one or more sets of two different n(z) samples use the form shear-shear=sample1-sample2 sample3-sample4 ....  Otherwise just use shear-shear=T to use the standard form.")


        #If no other spectra are specified, just do the shear-shear spectrum.
        if not any_spectra_option_found:
            print()
            print("No spectra requested in the parameter file so I will ")
            print("Assume you just want shear-shear, because it's the best one.")
            print()
            power3D = self.spectrumType.ShearShear.power_3d_type()
            shear_shear = self.spectrumType.ShearShear.value(self, power3D)
            self.req_spectra.append(shear_shear)
            self.req_power.append(power3D)
        elif not self.req_spectra:
            print()
            print("You switched off all the spectra in the parameter file")
            print("for project_2d.  I will go along with this and just do nothing,")
            print("but if you get a crash later this is probably why.")
            print()

    def load_distance_splines(self, block):
        #Extract some useful distance splines
        #have to copy these to get into C ordering (because we reverse them)
        z_distance = block[names.distances, 'z']
        a_distance = block[names.distances, 'a']
        chi_distance = block[names.distances, 'd_m']
        if z_distance[1]<z_distance[0]:
            z_distance = z_distance[::-1].copy()
            a_distance = a_distance[::-1].copy()
            chi_distance = chi_distance[::-1].copy()

        h0 = block[names.cosmological_parameters, "h0"]

        #convert Mpc to Mpc/h
        chi_distance *= h0
        
        if block.has_value(names.distances, 'CHISTAR'):
            self.chi_star = block[names.distances, 'CHISTAR'] * h0
        else:
            self.chi_star = None
        self.chi_max = chi_distance.max()
        self.a_of_chi = GSLSpline(chi_distance, a_distance)
        self.chi_of_z = GSLSpline(z_distance, chi_distance)

    #def load_rsd(self, block):
    """
    def load_matter_power(self, block, chi_of_z, k_growth=1.e-3):
        self.z_pk,self.k_pk,self.pk_lin = block.get_grid(names.matter_power_lin, 'z', 'k_h', 'p_k')
        growth_ind=np.where(self.k_pk>k_growth)[0][0]
        self.growth_array = np.sqrt(np.divide(self.pk_lin[:,growth_ind], self.pk_lin[0,growth_ind], 
                        out=np.zeros_like(self.pk_lin[:,growth_ind]), where=self.pk_lin[:,growth_ind]!=0.))
        self.chi_pk = chi_of_z(self.z_pk)

        _,self.k_pk,self.pk_nl = block.get_grid(names.matter_power_nl, 'z', 'k_h', 'p_k')

        log_growth = np.log(self.growth_array)
        log_a = np.log(1./(1+self.z_pk))
        log_a, log_growth = log_a[::-1], log_growth[::-1]
        logD_loga_spline = interp.InterpolatedUnivariateSpline(log_a, log_growth)
        self.f_of_z_pk = (logD_loga_spline.derivative()(log_a))[::-1]
        #print self.f_of_z_pk
    """

    def load_kernels(self, block):
        #During the setup we already decided what kernels (W(z) or N(z) splines)
        #we needed for the spectra we want to do.  Their names are stored in the
        #kernels_A and kernels_B dictionaries.
        for kernel_name, kernel_dict in self.kernels_A.items():
            #Check for any old kernels that should have been cleaned up by the
            #clean
            assert len(kernel_dict) == 0, "Internal cosmosis error: old cosmology not properly cleaned"
            #Load in the new kernel
            if self.verbose:
                print("Loading kernel {}".format(kernel_name))
            self.load_kernel(block, kernel_name, kernel_dict)

        #We are always using two kernels, which we call A and B
        #Often these will be the same groups of kernels, but not
        #always, so we may have to load them separately
        for kernel_name, kernel_dict in self.kernels_B.items():
            #Again, check for old kernels 
            assert len(kernel_dict) == 0, "Internal cosmosis error: old cosmology not properly cleaned"
            #Most of the time we are cross-correlating
            #samples and the kernel will be the same for A and B
            #In that case just refer to the one we already loaded
            if kernel_name in self.kernels_A:
            #Loop through all the loaded N(z), W(z)
                if self.verbose:
                    print("Already calculated {}".format(kernel_name))
                for i,kernel in self.kernels_A[kernel_name].items():
                    kernel_dict[i] = kernel
            else:
                #If not already cached we will have to load it freshly
                self.load_kernel(block, kernel_name, kernel_dict)

    def save_kernels(self, block, zmax):
        z = np.linspace(0.0, zmax, 1000)
        chi = self.chi_of_z(z)
        for kernel_name, kernel_dict in self.kernels_A.items() + self.kernels_B.items():
            section = "kernel_"+kernel_name
            block[section, "z"] = z
            block[section, "chi"] = chi
            for bin_i, kernel in kernel_dict.items():
                key = "bin_{}".format(bin_i)
                if not block.has_value(section, key):
                    nz = kernel(chi)
                    block[section,key] = nz

    def load_kernel(self, block, kernel_name, kernel_dict):
        #the name is of the form N_SAMPLENAME or W_SAMPLENAME, or K_SAMPLENAME
        kernel_type = kernel_name[0]
        sample_name = "nz_" + kernel_name[2:]
        if kernel_name[2:] == names.wl_number_density:
            sample_name = names.wl_number_density

        if kernel_type == "K":
            nbin = 1
        else:
            z = block[sample_name, 'z']
            nbin = block[sample_name, 'nbin']

        #Now load n(z) or W(z) for each bin in the range
        for i in range(nbin):
            if kernel_type=="N":
                kernel = limber.get_named_nchi_spline(block, sample_name, i+1, z, self.a_of_chi, self.chi_of_z)
            elif kernel_type=="W":
                kernel = limber.get_named_w_spline(block, sample_name, i+1, z, self.chi_max, self.a_of_chi)
            elif kernel_type=="K":
                if self.chi_star is None:
                    raise ValueError("Need to calculate chistar (comoving distance to last scattering) e.g. with camb to use CMB lensing.")
                curv_K = block.get_double(names.cosmological_parameters, "K", 0.0)
                kernel = limber.get_cmb_kappa_spline(self.chi_max, self.chi_star, self.a_of_chi, curv_K)
            else:
                raise ValueError("Unknown kernel type {0} ({1})".format(kernel_type, kernel_name))
            if kernel is None:
                raise ValueError("Could not load one of the W(z) or n(z) splines needed for limber integral (name={}, type={}, bin={})".format(kernel_name, kernel_type, i+1))
            kernel_dict[i] = kernel

    def load_one_power(self,  block, powerType):
        self.power[powerType], self.growth[powerType] = limber.load_power_growth_chi(
            block, self.chi_of_z, powerType.section_name, "k_h", "z", "p_k")

    def load_power(self, block):
        for powerType in self.req_power:
            self.load_one_power(block, powerType)

    def compute_spectra(self, block, spectrum):
        spectrum_name = spectrum.get_name()
        block[spectrum_name, 'ell'] = self.ell
        na, nb = spectrum.nbins()
        if spectrum.is_autocorrelation():
            block[spectrum_name, 'nbin'] = na
        block[spectrum_name, 'nbin_a'] = na
        block[spectrum_name, 'nbin_b'] = nb
        spectrum.prep_spectrum(block, self.chi_of_z, na, nbin2=nb, use_galaxy_power=self.use_galaxy_power)
        for i in range(na):
            #for auto-correlations C_ij = C_ji so we calculate only one of them,
            #but save both orderings to the block to account for different ordering
            #conventions.
            #for cross-correlations we must do both
            jmax = i+1 if spectrum.is_autocorrelation() else nb
            for j in range(jmax):
                c_ell = spectrum.compute(block, self.ell, i, j, 
                                        relative_tolerance=self.relative_tolerance, 
                                        absolute_tolerance=self.absolute_tolerance, 
                                        do_extended_limber=self.do_extended_limber,
                                        )
                block[spectrum_name, 'bin_{}_{}'.format(i+1,j+1)] = c_ell
                if self.get_kernel_peaks:
                    chi_peak, z_peak=spectrum.kernel_peak(block, i, j, self.a_of_chi)
                    block[spectrum_name, "chi_peak_{}_{}".format(i+1,j+1)] = chi_peak
                    block[spectrum_name, "z_peak_{}_{}".format(i+1,j+1)] = z_peak
                    block[spectrum_name, "arcmin_per_Mpch_{}_{}".format(i+1,j+1)] = 60*np.degrees(1/chi_peak)

    def clean(self):
        #need to manually delete power spectra we have loaded
        self.power.clear()
        self.growth.clear()

        #spectra know how to delete themselves, in gsl_wrappers.py
        for name,kernels in self.kernels_A.items():
            kernels.clear()
        for kernels in self.kernels_B.values():
            kernels.clear()
        self.outputs.clear()

    def execute(self, block):
        try:
            self.load_distance_splines(block)
            #self.load_matter_power(block, self.chi_of_z)
            try:
                self.load_kernels(block)
            except NullSplineError:
                sys.stderr.write("Failed to load one of the kernels (n(z) or W(z)) needed to compute 2D spectra\n")
                sys.stderr.write("Often this is because you are in a weird part of parameter space, but if it is \n")
                sys.stderr.write("consistent then you may wish to look into it in more detail. Set fata_errors=T to do so.\n")
                if self.fatal_errors:
                    raise
                return 1
            if self.save_kernel_zmax>0:
                self.save_kernels(block, self.save_kernel_zmax)
            self.load_power(block)
            for spectrum in self.req_spectra:
                if self.verbose:
                    print("Computing spectrum: {} -> {}".format(spectrum.__class__.__name__, spectrum.get_name()))
                self.compute_spectra(block, spectrum)
        finally:
            self.clean()
        return 0


def setup(options):
    return SpectrumCalculator(options)

def execute(block, config):
    return config.execute(block)
    
