"""

Errors in cosmic shear measurement can lead to a multiplicative factor
scaling the observed shear spectra.

This module scales the measured C_ell to account for that difference,
assuming model values of the multplicative factor m, either per bin or for all bins.

"""
from __future__ import print_function
from builtins import range
from cosmosis.datablock import names, option_section
import sys

warning_note_displayed = False


def setup(options):
    # This is an option - can set m_per_bin = T to get
    # a different m for each tomographic bin, or F to get
    # one global value
    m_per_bin = options.get_bool(option_section, "m_per_bin", True)
    cl_section = options.get_string(
        option_section, "cl_section", default=names.shear_cl)
    cross_section = options.get_string(
        option_section, "cross_section", default="galaxy_shear_cl")
    cmbcross_section=options.get_string(
        option_section, "cmbcross_section", default="shear_cmbkappa_cl")
    cal_section = options.get_string(
        option_section, "cal_section", default=names.shear_calibration_parameters)
    verbose = options.get_bool(option_section, "verbose", False)
    print()
    print("The shear_m_bias module will use calibration values from {} and look for ".format(cal_section))
    print("shear-shear spectra in {} and position-shear in {}".format(cl_section, cross_section))
    return m_per_bin, cl_section, cal_section, cross_section, cmbcross_section, verbose


def get_nbins(block, section):
    if block.has_value(section, "nbin_a"):
        n_a = block[section, "nbin_a"]
        n_b = block[section, "nbin_b"]
    else:
        n_a = block[section, "nbin"]
        n_b = n_a
    return n_a, n_b


def calibrate_section(block, section, m_a, m_b, verbose):
    n_a = len(m_a)
    n_b = len(m_b)
    for i in range(n_a):
        for j in range(n_b):

            # Get existing C_ell
            cl_name = "bin_{}_{}".format(i + 1, j + 1)
            if block.has_value(section, cl_name):
                if verbose:
                    print("Calibrating {} bin {} {} by (1+{}) * (1+{}) = {}".format(section, i + 1, j + 1, m_a[i], m_b[j], (1 + m_a[i]) * (1 + m_b[j])))
                block[section, cl_name] *= (1 + m_a[i]) * (1 + m_b[j])
            elif verbose:
                print("No {} bin {} {} to calibrate".format(section, i + 1, j + 1))


def calibrate_shear_shear(block, section, cal_section, m_per_bin, verbose):
    nbin_a, nbin_b = get_nbins(block, section)
    if m_per_bin:
        m = [block[cal_section, "m{}".format(i + 1)] for i in range(nbin_a)]
    else:
        m0 = block[cal_section, "m0"]
        m = [m0 for i in range(nbin_a)]
    calibrate_section(block, section, m, m, verbose)


def calibrate_position_shear(block, section, cal_section, m_per_bin, verbose):
    nbin_a, nbin_b = get_nbins(block, section)
    m_a = [0.0 for i in range(nbin_a)]
    if m_per_bin:
        m_b = [block[cal_section, "m{}".format(i + 1)] for i in range(nbin_b)]
    else:
        m0 = block[cal_section, "m0"]
        m_b = [m0 for i in range(nbin_b)]
    calibrate_section(block, section, m_a, m_b, verbose)

def calibrate_shear_cmbkappa(block, section, cal_section, m_per_bin, verbose):
    nbin_a, nbin_b = get_nbins(block, section)
    if (nbin_b != 1):
        raise IndexError("CMB kappa bins are set up incorrectly!")
    m_b = [0.0  for i in range(nbin_b)]
    if m_per_bin:
        m_a = [block[cal_section, "m{}".format(i+1)] for i in range(nbin_a)]
    else:
        m0 = block[cal_section, "m0"]
        m_a = [m0 for i in range(nbin_a)]
    calibrate_section(block, section, m_a, m_b, verbose)


def execute(block, config):
    m_per_bin, cl_section, cal_section, cross_section, cmbcross_section, verbose=config
    do_auto = block.has_section(cl_section)
    do_cross = block.has_section(cross_section)
    do_cmbcross = block.has_section(cmbcross_section)

    if do_auto:
        calibrate_shear_shear(
            block, cl_section, cal_section, m_per_bin, verbose)
    if do_cross:
        calibrate_position_shear(block, cross_section,
                                 cal_section, m_per_bin, verbose)
    if do_cmbcross:
        calibrate_shear_cmbkappa(block, cmbcross_section, cal_section, 
                                m_per_bin, verbose)

    if (not do_auto) and (not do_cross) and (not do_cmbcross):
        sys.stderr.write("ERROR: The shear bias calibration module could not find either a section {} or a {} to calibrate.\n".format(
            cl_section, cross_section))
        sys.stderr.write("The module therefore has nothing to do and considers this an error.  You may need to either change settings in the module or the precedng pipeline, or remove the module altogether\n")
        return 1

    global warning_note_displayed
    if not warning_note_displayed:
        warning_note_displayed = True
        if not do_auto:
            sys.stderr.write(
                "Note: No shear-shear section {} was found to calibrate. I did calibrate position-shear in {}.\n".format(cl_section, cross_section))
        elif not do_cross:
            sys.stderr.write(
                "Note: No position-shear section {} was found to calibrate. I did calibrate shear-shear in {}.\n".format(cross_section, cl_section))

    return 0
