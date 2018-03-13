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
from project_2d import SpectrumCalculator, Enum, Power3D, Spectrum, limber


class PPFLensingPower3D(Power3D):
    section = "ppf_modified_matter"
    source_specific = False

# Supported functions


class SpectrumTypePPF(Enum):
    class ShearShear(Spectrum):
        power_3d_type = PPFLensingPower3D
        kernels = "W W"
        autocorrelation = True
        name = names.shear_cl
        prefactor_power = 2


class SpectrumCalculatorPPF(SpectrumCalculator):
    spectrumType = SpectrumTypePPF

    def __init__(self, options):
        super(SpectrumCalculatorPPF, self).__init__(options)
        self.flatten_k = options.get_bool(option_section, "flatten_k", False)

    def load_power(self, block):
        for powerType in self.req_power:
            # Here we detail the ways we have to modify
            if isinstance(powerType, PPFLensingPower3D):
                self.power[powerType], self.growth[powerType] = load_power_growth_chi_ppf(
                    block, self.chi_of_z, names.matter_power_nl, "k_h", "z", "p_k", self.flatten_k)
            else:
                self.load_one_power(block, powerType)


def setup(options):
    # raise ValueError("Sorry - this module is broken.")
    return SpectrumCalculatorPPF(options)


def execute(block, config):
    return config.execute(block)



def load_power_growth_chi_ppf(block, chi_of_z, section, k_name, z_name, p_name, flatten_k, k_growth=1.e-3):
    z,k,p = block.get_grid(section, z_name, k_name, p_name)
    z1,k1,D = block.get_grid(names.post_friedmann_parameters, "z", "k_h", "D")

    # Make a growth function spline
    chi = chi_of_z(z)
    growth_spline = limber.growth_from_power(chi, k, p, k_growth)

    # Map the  D values onto the same k points as our P values
    # Make an interpolator for D(chi,logk)
    chi1 = chi_of_z(z1)
    D_spline = limber.GSLSpline2d(chi1, np.log(k1), D.T, spline_type=limber.BICUBIC)
    # k value our interpolator goes up to
    logk1_max = np.log(k1.max())

    # Array on same grid as our P(chi, k) that we will fill in 
    D_values = np.zeros_like(p)

    # For each k value 
    for i,ki in enumerate(k):
        logk = np.log(k[i])
        # If outside range then maybe flatten to max value
        if flatten_k and (logk>logk1_max):
            logk = logk1_max

        # Evaluate at all chi points
        D_i = D_spline(chi, logk)
        D_values[:,i] = D_i

    # Multiply by the D scaling factor.
    p *= D_values**2
    power_spline = limber.GSLSpline2d(chi, np.log(k), p.T, spline_type=limber.BICUBIC)

    return power_spline, growth_spline
