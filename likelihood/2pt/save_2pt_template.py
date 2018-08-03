"""
This module saves cosmosis output to a 2-pt file after interpolating it 
to specified ell values. It is useful for making simulations, but
is not yet fully supported so please use it with care and check 
the results carefully.

"""
from __future__ import print_function


from builtins import range
from builtins import object
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
import twopoint
from twopoint_cosmosis import type_table
import gaussian_covariance
twopt_like = __import__('2pt_like')  # don't start .py files with a number!
SpectrumInterp = twopt_like.SpectrumInterp


def setup(options):
    input_file = options[option_section, "input"]
    output_file = options[option_section, "output"]

    input_data = twopoint.TwoPointFile.from_fits(input_file)


    return (output_file, input_data)


def update_spectrum_from_block(block, spectrum):

    # The dictionary type_table stores the codes used in the FITS files for
    # the types of spectrum
    type_codes = (spectrum.type1.name, spectrum.type2.name)
    section, angle_name, bin_format = type_table[type_codes]

    splines = {}

    angle = spectrum.angle

    theory_angle = block[section, angle_name]

    # Convert to radians if needed
    if spectrum.angle_unit in ['arcmin', 'arcsec', 'deg']:
        print("Converting from {} to radians".format(spectrum.angle_unit))
        old_unit = twopoint.ANGULAR_UNITS[spectrum.angle_unit]
        new_unit = twopoint.ANGULAR_UNITS['rad']
        angle_with_units = angle * old_unit
        angle = angle_with_units.to(new_unit).value


    for i,(b1,b2,ang) in enumerate(zip(spectrum.bin1,spectrum.bin2, angle)):
        if (b1,b2) in splines:
            interp = splines[(b1,b2)]
        else:
            name_12 = bin_format.format(b1,b2)
            name_21 = bin_format.format(b2,b1)
            if block.has_value(section, name_12):
                theory = block[section, name_12]
            elif block.has_value(section, name_21) and (spectrum.type1==spectrum.type2):
                theory = block[section, name_21]
            else:
                raise ValueError("Missing section/name {} {}".format(section, name_12))
            splines[(b1,b2)] = interp = SpectrumInterp(theory_angle, theory)
        
        spectrum.value[i] = interp(ang)






def execute(block, config):
    output_file, input_data = config
    for spectrum in input_data.spectra:
        update_spectrum_from_block(block, spectrum)
    input_data.to_fits(output_file, clobber=True)

    return 0
