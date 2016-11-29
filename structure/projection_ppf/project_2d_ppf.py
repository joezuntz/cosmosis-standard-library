import sys
import os
import ctypes
import numpy as np
from cosmosis.datablock import names, option_section, BlockError

#The directory this file is in
dirname = os.path.split(__file__)[0]
#The directory for the vanilla projection code
vanilla_dir = os.path.join(dirname, os.path.pardir, 'projection')
sys.path.append(vanilla_dir)

#Now we have the path to import this
from project_2d  import SpectrumCalculator, Enum, PowerType, Spectrum, limber

class PowerTypePPF(Enum):
    ppf_modified_matter = "ppf_modified_matter"

#Supported functions
class SpectrumTypePPF(Enum):
    class ShearShear(Spectrum):
        power_spectrum = PowerTypePPF.ppf_modified_matter
        kernels = "W W"
        autocorrelation = True
        name = names.shear_cl
        prefactor_power = 2



class SpectrumCalculatorPPF(SpectrumCalculator):
    spectrumType = SpectrumTypePPF
    defaultSpectra=[spectrumType.ShearShear]

    def __init__(self, options):
        super(SpectrumCalculatorPPF,self).__init__(options)


    def load_ppf_power(self, block):
        #First load in the field D(k,z) which describes phi-psi
        MG_D = limber.load_power_chi(
            block, self.chi_of_z, names.post_friedmann_parameters, "k_h", "z", "D");


        #Make a function that will apply a scale to P(k,z) as 
        #we load it in.
        @limber.c_power_scaling_function
        def shear_shear_mg_scaling(k, z, P, args):
            D = limber.evaluate_power(MG_D, k, z)
            if D==0.0:
                return P
            else:
                return P * D**2
        return limber.load_power_chi_function(
                    block, self.chi_of_z, names.matter_power_nl, "k_h", "z", "p_k",
                    shear_shear_mg_scaling, None)

    def load_power(self, block):
        for powerType in self.req_power:
            #Here we detail the ways we have to modify 
            if powerType==PowerTypePPF.ppf_modified_matter:
                powerInterp = self.load_ppf_power(block)
            else:
                powerInterp = limber.load_power_chi(
                    block, self.chi_of_z, powerType.value, "k_h", "z", "p_k")
            self.power[powerType] = powerInterp


def setup(options):
    return SpectrumCalculatorPPF(options)

def execute(block, config):
    return config.execute(block)
    