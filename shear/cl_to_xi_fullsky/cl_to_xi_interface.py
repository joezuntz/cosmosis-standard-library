#coding: utf-8
#import cl_to_xi_full
from __future__ import print_function
from builtins import range
import numpy as np
from cosmosis.datablock import option_section, names as section_names
from cl_to_xi import save_xi_00_02, save_xi_22, arcmin_to_radians, SpectrumInterp, cl_to_xi_to_block
from legendre import get_legfactors_00, get_legfactors_02, get_legfactors_22, precomp_GpGm, apply_filter
from past.builtins import basestring

def setup(options):
    
    if options.has_value(option_section, "theta"):
        theta = options[option_section, 'theta']
        if np.isscalar(theta):
            theta = np.array([theta])
        theta = arcmin_to_radians(theta)
    else:
        n_theta = options.get_int(option_section, "n_theta", 20)
        theta_min = options.get_double(option_section, "theta_min", 1.) 
        theta_max = options.get_double(option_section, "theta_max", 300.)
        theta_min = arcmin_to_radians(theta_min)
        theta_max = arcmin_to_radians(theta_max)
        theta = np.logspace(np.log10(theta_min), np.log10(theta_max), n_theta)

    xi_type = options.get_string(option_section, 'xi_type')
    ell_max = options.get_int(option_section, "ell_max")
    #Filter the Cls at high ell to reduce ringing:
    high_l_filter = options.get_double(option_section, "high_l_filter", 0.75)

    cl_section = options.get_string(option_section, "input_section_name", "")
    output_section = options.get_string(option_section, "output_section_name", "")
    save_name = options.get_string(
        option_section, "save_name", "")

    # setup precompute functions and I/O sections

    if xi_type in ["22", "22+", "22-"]:
        precomp_func = get_legfactors_22
        if not cl_section:
            cl_section = "shear_cl"
        if output_section == "":
            output_section = (cl_section.replace("cl","xi_plus"), cl_section.replace("cl","xi_minus"))
        else:
            output_section = (output_section + "_plus", output_section + "_minus")
    elif xi_type == "00":
        precomp_func = get_legfactors_00
        if not cl_section:
            cl_section = "galaxy_cl"
        if output_section == "":
            output_section = cl_section.replace("cl", "xi")
    elif xi_type in ["02", "02+"]:
        precomp_func = get_legfactors_02
        if not cl_section:
            cl_section = "galaxy_shear_cl"
        if output_section == "":
            output_section = cl_section.replace("cl", "xi")
    else:
        print("xi_type should be '22', '22+' or '22-' for spin 2 autocorrelations e.g. xi+/-(theta),")
        print("00 for scalar autocorrelations e.g. w(theta)")
        print("or '02' or '02+' for spin 0 x spin 2 correlations e.g. gamma_t(theta)")
        raise ValueError()

    if save_name != "":
        add = "_%s"%save_name
        try:
            output_section += add
        except TypeError:
            output_section = ( output_section[0]+add, output_section[1]+add )

    legfacs = precomp_func(np.arange(ell_max + 1), theta)
    if high_l_filter>0:
        if isinstance(legfacs, tuple):
            legfacs = ( apply_filter( ell_max, high_l_filter, legfacs[0] ), 
            apply_filter( ell_max, high_l_filter, legfacs[1] ) )
        else:
            legfacs = apply_filter( ell_max, high_l_filter, legfacs )
    return theta, ell_max, legfacs, cl_section, output_section, save_name

def execute(block, config):

    thetas, ell_max, legfacs, cl_section, output_section, save_name = config

    ell = block[cl_section, "ell"]

    nbina, nbinb = block[cl_section, 'nbin_a'], block[cl_section, 'nbin_b']

    #Copy over info from 
    #block[output_section, "nbin_a"] = nbina
    #block[output_section, "nbin_b"] = nbinb
#
    #block[output_section, "theta"] = thetas
    #block[output_section, "save_name"] = save_name
    #block.put_metadata(output_section, "theta", "unit", "radians")

    if not isinstance(output_section, tuple):
        output_section = (output_section,)
        legfacs = (legfacs,)

    for i in range(1, nbina + 1):
        for j in range(1, nbinb + 1):
            name = 'bin_%d_%d' % (i, j)
            if block.has_value(cl_section, name):
                c_ell = block[cl_section, name]
            else:
                continue
            cl_interp = SpectrumInterp(ell, c_ell)
            cl_to_xi_to_block(block, output_section, name,
                          cl_interp, thetas, legfacs)

    if isinstance(output_section, basestring):
        output_section = (output_section,)
    for o in output_section:
        block[o, "nbin_a"] = nbina
        block[o, "nbin_b"] = nbinb
        block[o, "sample_a"] = block[cl_section, "sample_a"]
        block[o, "sample_b"] = block[cl_section, "sample_b"]
        block[o, "is_auto"] = block[cl_section, "is_auto"]
        block[o, "cl_section"] = cl_section
        block[o, "theta"] = thetas
        block[o, "sep_name"] = "theta"
        block[o, "save_name"] = save_name
        block.put_metadata(o, "theta", "unit", "radians")
    return 0

def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness.  The joy of python.
    return 0
