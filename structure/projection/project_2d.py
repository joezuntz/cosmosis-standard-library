#coding:utf-8
import os
import ctypes as ct
import numpy as np
import limber
from gsl_wrappers import GSLSpline
from cosmosis.datablock import names, option_section, BlockError
from enum34 import Enum
import re


class Spectrum(object):
    power_spectrum = "?"
    kernels = "??"
    autocorrelation = False
    name = "?"
    prefactor_power = np.nan

    def __init__(self, source, sample_a=names.wl_number_density, sample_b=names.wl_number_density):
        #caches of n(z), w(z), P(k,z), etc.
        self.source = source
        self.sample_a = sample_a
        self.sample_b = sample_b

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


# This is pretty cool.
# You can make an enumeration class which
# contains a list of possible options for something.
# But the options can each be anything - full python objects, and in this
# case they are classes of spectrum. So we can easily look these up by name,
# loop through them, etc.
class SpectrumType(Enum):

    class ShearShear(Spectrum):
        power_spectrum = "p_mm"
        kernels = "W W"
        autocorrelation = True
        name = names.shear_cl
        prefactor_power = 2

    class ShearIntrinsic(Spectrum):
        power_spectrum = "p_mi"
        kernels = "W N"
        autocorrelation = False
        name = names.shear_cl_gi
        prefactor_power = 1

    class IntrinsicIntrinsic(Spectrum):
        power_spectrum = "p_ii"
        kernels = "N N"
        autocorrelation = True
        name = names.shear_cl_ii
        prefactor_power = 0

    class PositionPosition(Spectrum):
        power_spectrum = "p_gg"
        kernels = "N N"
        autocorrelation = True
        name = "galaxy_cl"
        prefactor_power = 0

    class PositionMagnification(Spectrum):
        power_spectrum = "p_gm"
        kernels = "N W"
        autocorrelation = False
        name = "galaxy_magnification_cl"
        prefactor_power = 1    

    class MagnificationMagnification(Spectrum):
        power_spectrum = "p_mm"
        kernels = "W W"
        autocorrelation = True
        name = "magnification_cl"
        prefactor_power = 2

    class PositionShear(Spectrum):
        power_spectrum = "p_gm"
        kernels = "N W"
        autocorrelation = False
        name = "galaxy_shear_cl"
        prefactor_power = 1

    class PositionIntrinsic(Spectrum):
        power_spectrum = "p_gi"
        kernels = "N N"
        autocorrelation = False
        name = "galaxy_intrinsic_cl"
        prefactor_power = 0

    class MagnificationIntrinsic(Spectrum):
        power_spectrum = "p_mi"
        kernels = "W N"
        autocorrelation = False
        name = "magnification_intrinsic_cl"
        prefactor_power = 1
       
    class MagnificationShear(Spectrum):
        power_spectrum = "p_mm"
        kernels = "W W"
        autocorrelation = False
        name = "magnification_shear_cl"
        prefactor_power = 1




class SpectrumCalulcator(object):
    def __init__(self, options):
        #General options
        self.verbose = options.get_bool(option_section, "verbose", False)

        #Get the list of spectra that we want to compute.
        #The full list
        self.req_spectra = []
        for spectrum in SpectrumType:
            #By default we just do the shear-shear spectrum.
            #everything else is not done by default
            default = (spectrum==SpectrumType.ShearShear)
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
                #shear-shear = euclid-ska
                #where we would now search for nz_euclid and nz_ska
                values = value.split()
                for value in values:
                    try:
                        kernel_a, kernel_b = value.split('-',1)
                        kernel_a = kernel_a.strip()
                        kernel_b = kernel_b.strip()
                        self.req_spectra.append(spectrum.value(self, kernel_a, kernel_b))
                    except:
                        raise ValueError("To specify a P(k)->C_ell projection with one or more sets of two different n(z) samples use the form shear-shear=sample1-sample2 sample3-sample4 ....  Otherwise just use shear-shear=T to use the standard form.")



        print "Will project these spectra into 2D:"
        for spectrum in self.req_spectra:
            print "    - ", spectrum.name

        self.unbiased_galaxies = options.get_bool(option_section, 'unbiased_galaxies', default=False)

        if self.unbiased_galaxies:
            print "You have asked for 'unbiased_galaxies=T', so I will use the 3D matter power spectrum as the 3D galaxy power spectrum."

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
        
        self.chi_max = chi_distance.max()
        self.a_of_chi = GSLSpline(chi_distance, a_distance)
        self.chi_of_z = GSLSpline(z_distance, chi_distance)


    def load_kernels(self, block):
        #During the setup we already decided what kernels (W(z) or N(z) splines)
        #we needed for the spectra we want to do.  Their names are stored in the
        #kernels_A and kernels_B dictionaries.
        for kernel_name, kernel_dict in self.kernels_A.items():
            #Get rid of any N(z) from the previous point in the chain
            kernel_dict.clear()
            #Load in the new kernel
            if self.verbose:
                print "Loading kernel", kernel_name
            self.load_kernel(block, kernel_name, kernel_dict)

        #We are always using two kernels, which we call A and B
        #Often these will be the same groups of kernels, but not
        #always, so we may have to load them separately
        for kernel_name, kernel_dict in self.kernels_B.items():
            #Again, remove anything from the previous chain step
            kernel_dict.clear()
            #Most of the time we are cross-correlating
            #samples and the kernel will be the same for A and B
            #In that case just refer to the one we already loaded
            if kernel_name in self.kernels_A:
            #Loop through all the loaded N(z), W(z)
                if self.verbose:
                    print "Using cached ", kernel_name
                for i,kernel in self.kernels_A[kernel_name].items():
                    kernel_dict[i] = kernel
            else:
                #If not already cached we will have to load it freshly
                self.load_kernel(block, kernel_name, kernel_dict)



    def load_kernel(self, block, kernel_name, kernel_dict):
        #the name is of the form N_SAMPLENAME or W_SAMPLENAME
        kernel_type = kernel_name[0]
        sample_name = "nz_" + kernel_name[2:]
        if kernel_name[2:] == names.wl_number_density:
            sample_name = names.wl_number_density

        z = block[sample_name, 'z']
        nbin = block[sample_name, 'nbin']

        #Now load n(z) or W(z) for each bin in the range
        for i in xrange(nbin):
            if kernel_type=="N":
                kernel_dict[i] = limber.get_named_nchi_spline(block, sample_name, i+1, z, self.a_of_chi, self.chi_of_z)
            elif kernel_type=="W":
                kernel_dict[i] = limber.get_named_w_spline(block, sample_name, i+1, z, self.chi_max, self.a_of_chi)
            else:
                raise ValueError("Unknown kernel type {0} ({1})".format(kernel_type, kernel_name))


    def load_power(self, block):
        if 'p_mm' in self.req_power:
            self.power['p_mm'] = self.load_matter_power(block)

        if 'p_gg' in self.req_power:
            self.load_galaxy_power(block)

        if 'p_gm' in self.req_power:
            self.load_galaxy_matter_power(block)

        if 'p_mi' in self.req_power:
            self.power['p_mi'] = limber.load_power_chi(block, self.chi_of_z,
                                    names.intrinsic_alignment_parameters, "k_h", "z", "p_gi")
        if 'p_ii' in self.req_power:
            self.power['p_ii'] = limber.load_power_chi(block, self.chi_of_z,
                                    names.intrinsic_alignment_parameters, "k_h", "z", "p_ii")
        if 'p_gi' in self.req_power:
            self.load_galaxy_intrinsic_power(block)



    def load_matter_power(self, block):
        #If we have already loaded matter power just use that.
        #This one is a special case as we may use it elsewhere
        if 'p_mm' in self.power:
            p_mm = self.power['p_mm']
        #Otherwise load it.
        else:
            p_mm = limber.load_power_chi(block, self.chi_of_z,
                        names.matter_power_nl, "k_h", "z", "p_k")
        return p_mm


    def load_galaxy_power(self, block):
        #galaxy power. If we have been asked for unbiased galaxies we will
        #just use the matter power, either the one we have loaded
        if self.unbiased_galaxies:
            p_gg = self.load_matter_power(block)
        #we want the galaxy power directly
        else:
            #special case here for ease of use - give a more useful message
            #if galaxy power is not found
            if not block.has_section("galaxy_power"):
                raise RuntimeError("No galaxy_power section available. If you just want to use matter_power_nl then set unbiased_galaxies=T")
            #load the galaxy power spectrum
            p_gg = limber.load_power_chi(block, self.chi_of_z,
                'galaxy_power', "k_h", "z", "p_k")
        self.power['p_gg'] = p_gg

    def load_galaxy_matter_power(self, block):

        #galaxy power. If we have been asked for unbiased galaxies we will
        #just use the matter power, either the one we have loaded
        if self.unbiased_galaxies:
            p_gm = self.load_matter_power(block)
        else:
            #we want the galaxy power directly
            #special case here for ease of use - give a more useful message
            #if galaxy power is not found
            if not block.has_section("galaxy_matter_power"):
                raise RuntimeError("No galaxy_power section available. If you just want to use matter_power_nl then set unbiased_galaxies=T")
            #load the galaxy power spectrum
            p_gm = limber.load_power_chi(block, self.chi_of_z,
                'galaxy_matter_power', "k_h", "z", "p_k")
        self.power['p_gm'] = p_gm

    def load_galaxy_intrinsic_power(self, block):
        #galaxy power. If we have been asked for unbiased galaxies we will
        #just use the matter power, either the one we have loaded
        if self.unbiased_galaxies:
            if 'p_mi' in self.power:
                p_gi = self.power['p_mi']
            else:
                p_gi = limber.load_power_chi(block, self.chi_of_z,
                    names.intrinsic_alignment_parameters, "k_h", "z", "p_mi")
        else:
            #load the galaxy power spectrum
            p_gi = limber.load_power_chi(block, self.chi_of_z,
                'names.intrinsic_alignment_parameters', "k_h", "z", "p_pos_i")
        self.power['p_gi'] = p_gi


    def compute_spectra(self, block, spectrum):
        block[spectrum.name, 'ell'] = self.ell
        na, nb = spectrum.nbins()
        for i in xrange(na):
            #for auto-correlations C_ij = C_ji so we only do one of them.
            #for cross-correlations we must do both
            jmax = i+1 if spectrum.is_autocorrelation() else nb
            for j in xrange(jmax):
                c_ell = spectrum.compute(block, self.ell, i, j)
                self.outputs[spectrum.name+"_{}_{}".format(i,j)] = c_ell
                block[spectrum.name, 'bin_{}_{}'.format(i+1,j+1)] = c_ell(self.ell)
            

    def clean(self):
        self.power.clear()
        for kernels in self.kernels_A.values():
            kernels.clear()
        for kernels in self.kernels_B.values():
            kernels.clear()
        self.outputs.clear()

    def execute(self, block):
        self.load_distance_splines(block)
        self.load_kernels(block)
        self.load_power(block)
        for spectrum in self.req_spectra:
            print "Computing spectrum:", spectrum.__class__.__name__
            self.compute_spectra(block, spectrum)
        self.clean()
        return 0




def setup(options):
    return SpectrumCalulcator(options)

def execute(block, config):
    return config.execute(block)
    