#coding: utf-8
#import cl_to_xi_full
from __future__ import print_function
from builtins import range
import numpy as np
from astropy.units import Quantity
from cosmosis.datablock import option_section, names as section_names


def setup(options):

    output_units = options.get_string(option_section, 'output_units')
    section = options.get_string(option_section, 'section_name', '')

    if output_units not in ['rad', 'deg', 'arcmin', 'arcsec']:
        raise ValueError('output units need to be one of {}'.format(['rad', 'deg', 'arcmin', 'arcsec']))

    return output_units, section


def execute(block, config):

    output_units, section = config
    
    input_units = block.get_metadata(section, "theta", "unit")
    theta_in = block[section, "theta"]
    if input_units in ['radian', 'radians']:
        input_units = 'rad'
    
    theta_in_astropy = Quantity(theta_in, input_units)
    theta_out_astropy = theta_in_astropy.to(output_units)

    block.replace_double_array_1d(section, "theta", theta_out_astropy.value)
    block.replace_metadata(section, "theta", "unit", output_units)
    
    return 0

def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness.  The joy of python.
    return 0
