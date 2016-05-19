#coding:utf-8
import os
import ctypes as ct
import numpy as np
import limber
from gsl_wrappers import GSLSpline
from cosmosis.datablock import names, option_section, BlockError
from enum34 import Enum
import re




class PowerType(Enum):
    matter = "matter_power_nl" #special case
    galaxy = "galaxy_power"
    intrinsic = "intrinsic_power"
    matter_galaxy = "matter_galaxy_power"
    matter_intrinsic = "matter_intrinsic_power"
    galaxy_intrinsic = "galaxy_intrinsic_power"




class Spectrum(object):
    power_spectrum = "?"
    kernels = "??"
    autocorrelation = False
    name = "?"
    prefactor_power = np.nan

    def __init__(self, source, sample_a=names.wl_number_density, sample_b=names.wl_number_density, save_name=""):
        #caches of n(z), w(z), P(k,z), etc.
        self.source = source
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


    def prefactor(self, block):
        if self.prefactor_power == 0:
            return 1.0
        c_kms = 299792.4580
        omega_m = block[names.cosmological_parameters, "omega_m"]
        shear_scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m
        return shear_scaling**self.prefactor_power
    
    def magnification_prefactor(self, block, bin):
        alpha = block[names.galaxy_luminosity_function, "alpha_binned"]
        return alpha[bin]


    def compute(self, block, ell, bin1, bin2):
        #Get the required kernels
        #maybe get from rescaled version of existing spectra
        P = self.source.power[self.power_spectrum]
        K1 = self.source.kernels_A[self.kernels[ 0]+"_"+self.sample_a][bin1]
        K2 = self.source.kernels_B[self.kernels[-1]+"_"+self.sample_b][bin2]

        #The spline is done in log space if the spectrum never goes
        #negative, otherwise linear space.
        xlog = ylog = self.autocorrelation
        #Generate a spline
        c_ell = limber.limber(K1, K2, P, xlog, ylog, ell, self.prefactor(block))
        return c_ell

    def kernel_peak(self, block, bin1, bin2):
        K1 = self.source.kernels_A[self.kernels[ 0]+"_"+self.sample_a][bin1]
        K2 = self.source.kernels_B[self.kernels[-1]+"_"+self.sample_b][bin2]
        return limber.get_kernel_peak(K1, K2)

# This is pretty cool.
# You can make an enumeration class which
# contains a list of possible options for something.
# But the options can each be anything - full python objects, and in this
# case they are classes of spectrum. So we can easily look these up by name,
# loop through them, etc.
class SpectrumType(Enum):

    class ShearShear(Spectrum):
        power_spectrum = PowerType.matter
        kernels = "W W"
        autocorrelation = True
        name = names.shear_cl
        prefactor_power = 2

    class ShearIntrinsic(Spectrum):
        power_spectrum = PowerType.matter_intrinsic
        kernels = "W N"
        autocorrelation = False
        name = names.shear_cl_gi
        prefactor_power = 1

    class IntrinsicIntrinsic(Spectrum):
        power_spectrum = PowerType.intrinsic
        kernels = "N N"
        autocorrelation = True
        name = names.shear_cl_ii
        prefactor_power = 0

    class PositionPosition(Spectrum):
        power_spectrum = PowerType.galaxy
        kernels = "N N"
        autocorrelation = True
        name = "galaxy_cl"
        prefactor_power = 0

    class MagnificationPosition(Spectrum):
        power_spectrum = PowerType.matter_galaxy
        kernels = "W N"
        autocorrelation = False
        name = "magnification_galaxy_cl"
        prefactor_power = 1    

    class MagnificationMagnification(Spectrum):
        power_spectrum = PowerType.matter
        kernels = "W W"
        autocorrelation = True
        name = "magnification_cl"
        prefactor_power = 2

    class PositionShear(Spectrum):
        power_spectrum = PowerType.matter_galaxy
        kernels = "N W"
        autocorrelation = False
        name = "galaxy_shear_cl"
        prefactor_power = 1

    """
    class ShearPosition(Spectrum):
        power_spectrum = PowerType.matter_galaxy
        kernels = "W N"
        autocorrelation = False
        name = "shear_galaxy_cl"
        prefactor_power = 1
    """

    class PositionIntrinsic(Spectrum):
        power_spectrum = PowerType.galaxy_intrinsic
        kernels = "N N"
        autocorrelation = False
        name = "galaxy_intrinsic_cl"
        prefactor_power = 0

    class MagnificationIntrinsic(Spectrum):
        power_spectrum = PowerType.matter_intrinsic
        kernels = "W N"
        autocorrelation = False
        name = "magnification_intrinsic_cl"
        prefactor_power = 1
       
    class MagnificationShear(Spectrum):
        power_spectrum = PowerType.matter
        kernels = "W W"
        autocorrelation = False
        name = "magnification_shear_cl"
        prefactor_power = 1

    class ShearCmbkappa(Spectrum):
        power_spectrum = PowerType.matter
        kernels = "W K"
        autocorrelation = False
        name = "shear_cmbkappa_cl"
        prefactor_power = 2

    class CmbkappaCmbkappa(Spectrum):
        power_spectrum = PowerType.matter
        kernels = "K K"
        autocorrelation = True
        name = "cmbkappa_cl"
        prefactor_power = 2

    class IntrinsicCmbkappa(Spectrum):
        power_spectrum = PowerType.matter_intrinsic
        kernels = "N K"
        autocorrelation = False
        name = "intrinsic_cmbkappa_cl"
        prefactor_power = 1

    class PositionCmbkappa(Spectrum):
        power_spectrum = PowerType.matter_galaxy
        kernels = "N K"
        autocorrelation = False
        name = "galaxy_cmbkappa_cl"
        prefactor_power = 1





class SpectrumCalulcator(object):
    # It is useful to put this here so we can subclass to add new spectrum
    # types, for example ones done with modified gravity changes.
    spectrumType = SpectrumType
    def __init__(self, options):
        #General options
        self.verbose = options.get_bool(option_section, "verbose", False)
        self.get_kernel_peaks = options.get_bool(option_section, "get_kernel_peaks", False)

        #Get the list of spectra that we want to compute.
        #The full list
        self.req_spectra = []
        for spectrum in self.spectrumType:
            #By default we just do the shear-shear spectrum.
            #everything else is not done by default
            default = (spectrum==self.spectrumType.ShearShear)
            name = spectrum.value.option_name()

            try:
                #first look for a string, e.g. redmagic-redmagic or similar (name of samples)
                value = options.get_string(option_section, name)
            except BlockError:
                #or try a bool
                value = options.get_bool(option_section, name, default=default)
            if isinstance(value, bool) and value:
                self.req_spectra.append(spectrum.value(self))
                #Now we just look for the wl_number_density section.
                #though we may change this
            elif isinstance(value, bool):
                pass
            else:
                #now we are looking for things of the form
                #shear-shear = euclid-ska:
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
                        if ":" in kernel_b:
                            kernel_b, save_name=kernel_b.split(":",1)
                        else:
                            save_name = ""
                        kernel_a = kernel_a.strip()
                        kernel_b = kernel_b.strip()
                        #The self in the line below is not a mistake - the source objects
                        #for the spectrum class is the SpectrumCalculator itself
                        self.req_spectra.append(spectrum.value(self, kernel_a, kernel_b, save_name))
                    except:
                        raise ValueError("To specify a P(k)->C_ell projection with one or more sets of two different n(z) samples use the form shear-shear=sample1-sample2 sample3-sample4 ....  Otherwise just use shear-shear=T to use the standard form.")



        print "Will project these spectra into 2D:"
        for spectrum in self.req_spectra:
            print "    - ", spectrum.get_name()


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

        self.req_power = set()
        for spectrum in self.req_spectra:
            self.req_power.add(spectrum.power_spectrum)


        self.power = {}
        self.outputs = {}

        #And the req ell ranges.
        #We use log-spaced output
        ell_min = options.get_double(option_section, "ell_min")
        ell_max = options.get_double(option_section, "ell_max")
        n_ell = options.get_int(option_section, "n_ell")
        self.ell = np.logspace(np.log10(ell_min), np.log10(ell_max), n_ell)


    def load_distance_splines(self, block):
        #Extract some useful distance splines
        #have to copy these to get into C ordering (because we reverse them)
        z_distance = block[names.distances, 'z'][::-1].copy()
        a_distance = block[names.distances, 'a'][::-1].copy()
        chi_distance = block[names.distances, 'd_m'][::-1].copy()
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
                print "Loading kernel", kernel_name
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
                    print "Already calculated ", kernel_name
                for i,kernel in self.kernels_A[kernel_name].items():
                    kernel_dict[i] = kernel
            else:
                #If not already cached we will have to load it freshly
                self.load_kernel(block, kernel_name, kernel_dict)



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
        for i in xrange(nbin):
            if kernel_type=="N":
                kernel_dict[i] = limber.get_named_nchi_spline(block, sample_name, i+1, z, self.a_of_chi, self.chi_of_z)
            elif kernel_type=="W":
                kernel_dict[i] = limber.get_named_w_spline(block, sample_name, i+1, z, self.chi_max, self.a_of_chi)
            elif kernel_type=="K":
                if self.chi_star is None:
                    raise ValueError("Need to calculate chistar (comoving distance to last scattering) e.g. with camb to use CMB lensing.")
                kernel_dict[i] = limber.get_cmb_kappa_spline(self.chi_max, self.chi_star, self.a_of_chi)
            else:
                raise ValueError("Unknown kernel type {0} ({1})".format(kernel_type, kernel_name))


    def load_power(self, block):
        for powerType in self.req_power:
            self.power[powerType] = limber.load_power_chi(
                block, self.chi_of_z, powerType.value, "k_h", "z", "p_k")


    def compute_spectra(self, block, spectrum):
        spectrum_name = spectrum.get_name()
        block[spectrum_name, 'ell'] = self.ell
        na, nb = spectrum.nbins()
        if spectrum.is_autocorrelation():
            block[spectrum_name, 'nbin'] = na
        block[spectrum_name, 'nbin_a'] = na
        block[spectrum_name, 'nbin_b'] = nb
        for i in xrange(na):
            #for auto-correlations C_ij = C_ji so we calculate only one of them,
            #but save both orderings to the block to account for different ordering
            #conventions.
            #for cross-correlations we must do both
            jmax = i+1 if spectrum.is_autocorrelation() else nb
            for j in xrange(jmax):
                c_ell = spectrum.compute(block, self.ell, i, j)
                self.outputs[spectrum_name+"_{}_{}".format(i,j)] = c_ell
                block[spectrum_name, 'bin_{}_{}'.format(i+1,j+1)] = c_ell(self.ell)
                if self.get_kernel_peaks:
                    chi_peak=spectrum.kernel_peak(block, i, j)
                    block[spectrum_name, "chi_peak_{}_{}".format(i+1,j+1)] = chi_peak
                    block[spectrum_name, "arcmin_per_Mpch_{}_{}".format(i+1,j+1)] = 60*np.degrees(1/chi_peak)

    def clean(self):
        #need to manually delete power spectra we have loaded
        for p in self.power.values():
            limber.free_power(p)
        self.power.clear()

        #spectra know how to delete themselves, in gsl_wrappers.py
        for name,kernels in self.kernels_A.items():
            kernels.clear()
        for kernels in self.kernels_B.values():
            kernels.clear()
        self.outputs.clear()

    def execute(self, block):
        try:
            self.load_distance_splines(block)
            self.load_kernels(block)
            self.load_power(block)
            for spectrum in self.req_spectra:
                if self.verbose:
                    print "Computing spectrum: {} -> {}".format(spectrum.__class__.__name__, spectrum.get_name())
                self.compute_spectra(block, spectrum)
        finally:
            self.clean()
        return 0




def setup(options):
    return SpectrumCalulcator(options)

def execute(block, config):
    return config.execute(block)
    
