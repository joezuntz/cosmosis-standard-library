import sys
import os
import ctypes
import numpy as np
from cosmosis.datablock import names, option_section, BlockError

# The directory this file is in
dirname = os.path.split(__file__)[0]
# The directory for the vanilla projection code
vanilla_dir = os.path.join(dirname, os.path.pardir, 'projection')
sys.path.append(vanilla_dir)

# Now we have the path to import this
from project_2d import SpectrumCalculator, Enum, PowerType, Spectrum, limber


class PowerTypePPF(Enum):
    ppf_modified_matter = "ppf_modified_matter"

# Supported functions


class SpectrumTypePPF(Enum):
    class ShearShear(Spectrum):
        power_spectrum = PowerTypePPF.ppf_modified_matter
        kernels = "W W"
        autocorrelation = True
        name = names.shear_cl
        prefactor_power = 2


class SpectrumCalculatorPPF(SpectrumCalculator):
    spectrumType = SpectrumTypePPF
    defaultSpectra = [spectrumType.ShearShear]

    def __init__(self, options):
        super(SpectrumCalculatorPPF, self).__init__(options)

    def load_ppf_power(self, block):
        for powerType in self.req_power:
            self.power[powerType], self.growth[powerType] = load_power_growth_chi_ppf(
                block, self.chi_of_z, names.matter_power_nl, "k_h", "z", "p_k")

    def load_power(self, block):
        for powerType in self.req_power:
            # Here we detail the ways we have to modify
            if powerType == PowerTypePPF.ppf_modified_matter:
                powerInterp = self.load_ppf_power(block)
            else:
                powerInterp = super(SpectrumCalculatorPPF,self).load_power()
            self.power[powerType] = powerInterp


def setup(options):
    raise ValueError("Sorry - this module is broken.")
    return SpectrumCalculatorPPF(options)


def execute(block, config):
    return config.execute(block)



def load_power_growth_chi_ppf(block, chi_of_z, section, k_name, z_name, p_name, k_growth=1.e-3):
    z,k,p = block.get_grid(section, z_name, k_name, p_name)
    z1,k1,D = block.get_grid(names.post_friedmann_parameters, "z", "k_h", "D")

    chi = chi_of_z(z)
    growth_spline = limber.growth_from_power(chi, k, p, k_growth)

    # Map the  D values onto the same k points as our P values
    # Doesn't currently work because we need to extrapolate D(k) in general.
    chi1 = chi_of_z(z1)
    D_spline = limber.GSLSpline2d(chi1, np.log(k1), D.T, spline_type=limber.BICUBIC)
    D_values = np.zeros_like(p)
    for i in xrange(len(z)):
        D_i = D_spline(chi, np.log(k[i]))
        D_values[:,i] = D_i


    # Multiply by the D scaling factor.
    p *= D_values
    power_spline = limber.GSLSpline2d(chi, np.log(k), p.T, spline_type=limber.BICUBIC)
    return power_spline, growth_spline
