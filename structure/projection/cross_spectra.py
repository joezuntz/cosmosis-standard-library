#This is based on project_2d.py, but generalised for cross-spectra e.g. when correlating
#two populations with different n(z)s.

#coding:utf-8
import os
import ctypes as ct
import numpy as np
import limber_named as limber
from gsl_wrappers import GSLSpline
from cosmosis.datablock import names, option_section
from enum34 import Enum

class Spectrum(object):
    power_spectrum = "?"
    kernel_1 = "?"
    kernel_2 = "?"
    strictly_positive = False
    name = "?"
    prefactor_power = 1

    def __init__(self, source):
        #caches of n(z), w(z), P(k,z), etc.
        self.source = source

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
        K1 = self.source.kernels_A[self.kernels[0]][bin1]
        K2 = self.source.kernels_B[self.kernels[-1]][bin2]
        #The spline is done in log space if the spectrum never goes
        #negative, otherwise linear space.
        xlog = ylog = self.strictly_positive
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
        strictly_positive = True
        name = names.shear_cl
        prefactor_power = 2

    class ShearIntrinsic(Spectrum):
        power_spectrum = "p_mi"
        kernels = "W N"
        strictly_positive = False
        name = names.shear_cl_gi
        prefactor_power = 1

    class IntrinsicIntrinsic(Spectrum):
        power_spectrum = "p_ii"
        kernels = "N N"
        strictly_positive = True
        name = names.shear_cl_ii
        prefactor_power = 0

    class PositionPosition(Spectrum):
        power_spectrum = "p_gg"
        kernels = "N N"
        strictly_positive = True
        name = "galaxy_cl"
        prefactor_power = 0

    class PositionMagnification(Spectrum):
        power_spectrum = "p_gm"
        kernels = "N W"
        strictly_positive = False
        name = "galaxy_magnification_cl"
        prefactor_power = 1    

    class MagnificationMagnification(Spectrum):
        power_spectrum = "p_mm"
        kernels = "W W"
        strictly_positive = True
        name = "magnification_cl"
        prefactor_power = 2

    class PositionShear(Spectrum):
        power_spectrum = "p_gm"
        kernels = "N W"
        strictly_positive = False
        name = "galaxy_shear_cl"
        prefactor_power = 1

    class PositionIntrinsic(Spectrum):
        power_spectrum = "p_gi"
        kernels = "N N"
        strictly_positive = False
        name = "galaxy_intrinsic_cl"
        prefactor_power = 0

    class MagnificationIntrinsic(Spectrum):
        power_spectrum = "p_mi"
        kernels = "W N"
        strictly_positive = False
        name = "magnification_intrinsic_cl"
        prefactor_power = 1
       
    class MagnificationShear(Spectrum):
        power_spectrum = "p_mm"
        kernels = "W W"
        strictly_positive = False
        name = "magnification_shear_cl"
        prefactor_power = 1

class Kernels(object):
    def __init__(self, keys=["N","W"]):
        self.N={}
        self.W={}

class SpectrumCalulcator(object):
    def __init__(self, options):
        #General options
        self.verbose = options.get_bool(option_section, "verbose", False)

        self.nofz_sec_A=options.get_string(option_section, "nofz_section_A", names.wl_number_density)
        self.nofz_sec_B=options.get_string(option_section, "nofz_section_B", names.wl_number_density)



        #Get the list of spectra that we want to compute.
        #The full list
        self.req_spectra = []
        for spectrum in SpectrumType:
            #By default we just do the shear-shear spectrum.
            #everything else is not done by default
            default = (spectrum==SpectrumType.ShearShear)
            if options.get_bool(option_section, spectrum.value.name, default=default):
                self.req_spectra.append(spectrum.value(self))

        #Read section_prefix - if more than one cross-spectra in pipeline 
        #e.g. source-source clustering and source-lens clustering, these need
        #to be saved in different sections, so use this prefix
        self.section_prefix=options.get_string(option_section, "section_prefix", "None")
        if self.section_prefix!="None":
            for spectrum in self.req_spectra:
                spectrum.name = self.section_prefix+'_'+spectrum.name

        print "Will project these spectra into 2D:"
        for spectrum in self.req_spectra:
            print "    - ", spectrum.name

        self.unbiased_galaxies = options.get_bool(option_section, 'unbiased_galaxies', default=False)

        if self.unbiased_galaxies:
            print "You have asked for 'unbiased_galaxies=T', so I will use the 3D matter power spectrum as the 3D galaxy power spectrum."

        #Decide which kernels we will need to save
        self.req_kernels_A, self.req_kernels_B = set(),set()
        for spectrum in self.req_spectra:
            self.req_kernels_A.add(spectrum.kernels[0])
            self.req_kernels_B.add(spectrum.kernels[-1])

        self.req_power = set()
        for spectrum in self.req_spectra:
            self.req_power.add(spectrum.power_spectrum)


        #Modify this to do cross-experiment or sample 
        #spectra
        self.kernels_A = {"N":{}, "W":{}}
        self.kernels_B = {"N":{}, "W":{}}
        self.power = {}
        self.outputs = {}

        #For convenience in loading kernels, make some lists
        self.kernels = [self.kernels_A, self.kernels_B]
        self.req_kernels = [self.req_kernels_A, self.req_kernels_B]
        self.nofz_secs = [self.nofz_sec_A, self.nofz_sec_B]

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
        for kernels, req_kernels, nofz_sec, nbin in zip(self.kernels,self.req_kernels,self.nofz_secs,self.nbins):
            if 'N' in req_kernels:
                self.load_nchi_splines(block, kernels, nbin, nofz_sec)
            if 'W' in req_kernels:
                self.load_w_splines(block, kernels, nbin, nofz_sec)

    def load_nchi_splines(self, block, kernels, nbin, nofz_sec):
        kernels['N'].clear()
        z = block[nofz_sec, 'z']
        for i in xrange(nbin):
            kernels['N'][i] = limber.get_nchi_spline(block, nofz_sec, i+1, z, self.a_of_chi, self.chi_of_z)

    def load_w_splines(self, block, kernels, nbin, nofz_sec):
        kernels['W'].clear()
        z = block[nofz_sec, 'z']
        for i in xrange(nbin):
            kernels['W'][i] = limber.get_w_spline(block, nofz_sec, i+1, z, self.chi_max, self.a_of_chi)

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
        #For auto-correlations we don't need e.g. Cl_1_2 and Cl_2_1 - they're the same
        #For cross-correlations we do
        #We could add another option to the spectrum class for this...however the 
        #strictly_position option actually does this already - only auto-correlations are
        #always positive.
        for i in xrange(self.nbin_A):
            for j in xrange(self.nbin_B):
                c_ell = spectrum.compute(block, self.ell, i, j)
                self.outputs[spectrum.name+"_{}_{}".format(i,j)] = c_ell
                #if self.section_prefix is not None:
                #    output_section=self.section_prefix+"_"+spectrum.name
                #else:
                output_section=spectrum.name
                block[output_section, 'bin_{}_{}'.format(i+1,j+1)] = c_ell(self.ell)
            
    def clean(self):
        self.power.clear()
        for kernels in self.kernels_A.values():
            kernels.clear()
        self.outputs.clear()

    def execute(self, block):
        self.nbin_A = block[self.nofz_sec_A, 'nbin']
        self.nbin_B = block[self.nofz_sec_B, 'nbin']
        self.nbins = [self.nbin_A,self.nbin_B]
        self.load_distance_splines(block)
        self.load_kernels(block)
        self.load_power(block)
        for spectrum in self.req_spectra:
            print "Computing spectrum:", spectrum.__class__.__name__
            self.compute_spectra(block, spectrum)
            block[spectrum.name, "nbin_A"]=self.nbin_A
            block[spectrum.name, "nbin_B"]=self.nbin_B
        self.clean()
        return 0

def setup(options):

    return SpectrumCalulcator(options)

def execute(block, config):
    return config.execute(block)
    
