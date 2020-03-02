#coding: utf-8
#import cl_to_xi_full
from __future__ import print_function
from builtins import range
import numpy as np
from cosmosis.datablock import option_section, names as section_names
from cl_to_xi import save_xi_00_02, save_xi_22, arcmin_to_radians, SpectrumInterp, cl_to_xi_to_block
from legendre import get_legfactors_00, get_legfactors_02, get_legfactors_22, precomp_GpGm, apply_filter, get_legfactors_02_binav, get_legfactors_00_binav, get_legfactors_22_binav
from past.builtins import basestring
import sys
import os
dirname = os.path.split(__file__)[0]
twopoint_path = os.path.join(dirname,"..","..","likelihood","2pt")
sys.path.append(twopoint_path)
import twopoint

def readtheta(filename, xi_type_2pt, theta_type = 'centers', desired_units = 'arcmin'):
    '''
    Short function to read in theta values from a specified fits file.
    Desired angle units in 'rad', 'arcmin', 'arcsec', 'deg
    '''
    
    T = twopoint.TwoPointFile.from_fits(filename)
    xi = T.get_spectrum(xi_type_2pt)
    # make sure the units are in arcmin (or whatever else you want)
    xi.convert_angular_units(desired_units)
    # get the theta values for a single bin pair.
    # This method currently assumes that the same binning is used for all bin pairs (in a given corr funct).
    # scale cuts are applied later and can be different for each bin pair. 
    pairs = xi.get_bin_pairs()
    indices = xi.get_pair_mask(pairs[0][0], pairs[0][1]) #use the first appearing pair, since all should be the same.
    if theta_type == 'centers':
        thetas = xi.angle[indices]
        return thetas
    elif theta_type == 'edges':
        theta_min = xi.angle_min[indices]  # array of min values
        theta_max = xi.angle_max[indices] # array of max values
        # We currently assume adjacent and non-overlapping bins.
        np.testing.assert_allclose(theta_min[1:], theta_max[:-1], rtol=1e-03,
        err_msg='theta bin min and max values do not match.', verbose=True)
        theta_bins = np.append(theta_min, theta_max[-1])
        return theta_bins
    else:
        #if some other option is passed, eg "all", return both.
        thetas = xi.angle[indices]
        theta_min = xi.angle_min[indices]  # array of min values
        theta_max = xi.angle_max[indices] # array of max values
        np.testing.assert_allclose(theta_min[1:], theta_max[:-1], rtol=1e-03,
        err_msg='theta bin min and max values do not match.', verbose=True)
        theta_bins = np.append(theta_min, theta_max[-1])
        return thetas, theta_bins

def setup(options):
    
    xi_type_conversion = {'xip':['22', '22+'], 'xim':['22-'], 'gammat':['02', '02+'], 'wtheta':['00']}
    xi_type = options.get_string(option_section, 'xi_type')
    xi_type_2pt = None
    for key, value in xi_type_conversion.items():
        if xi_type in value:
            xi_type_2pt = key
            # check that this syntax works.
            break
    if xi_type_2pt is None:
        raise ValueError('xi_type %s not associated with a fits key'%xi_type)
    ell_max = options.get_int(option_section, "ell_max")
    bin_avg = options.get_bool(option_section, "bin_avg", False)
    theta_from_dat = options.get_bool(option_section, "theta_from_dat", False)
    two_thirds_midpoint = options.get_bool(option_section, "two_thirds_midpoint", True)
    #Filter the Cls at high ell to reduce ringing:
    high_l_filter = options.get_double(option_section, "high_l_filter", 0.75)

    cl_section = options.get_string(option_section, "input_section_name", "")
    output_section = options.get_string(option_section, "output_section_name", "")
    save_name = options.get_string(
        option_section, "save_name", "")
    theta_edges = None
    print('theta_from_dat',theta_from_dat)
    print('bin_avg',bin_avg)
    if theta_from_dat:
        print("*** Reading in theta values from datafile ***")
        if options.has_value(option_section, "filename"):
            filename = options.get_string(option_section, 'filename')
            print('filename = ',filename)
        else:
            raise ValueError('No filename specified for theta values')
    if bin_avg:
        print("*** Using bin averaged Legendre coefficients ***")
        if theta_from_dat:
            theta_edges = readtheta(filename, xi_type_2pt, theta_type = 'edges', desired_units = 'rad')
        else:
            if options.has_value(option_section, "theta"):
                print('WARNING: Specify `theta_edges` instead of `theta` when using bin averaging. `theta` is ignored')
                #raise ValueError()#maybe we don't want to actually raise an error, since we may be inheriting a default params.ini file.
            if options.has_value(option_section, "n_theta"):
                print('WARNING: Specify `n_theta_bins` instead of `n_theta` when using bin averaging. `n_theta` is ignored')
               # raise ValueError()
            if options.has_value(option_section, "theta_edges"):
                print('*** Note: Using `theta_edges` values specified in ini file ***')
                theta_edges = options[option_section, 'theta_edges']
                if np.isscalar(theta_edges):
                    theta_edges = np.array([theta_edges])
                theta_edges = arcmin_to_radians(theta_edges)
            else:
                print('** NOTE that theta_min and theta_max refer to the min and max of the first and last angular bin, respectively **')
                n_theta_bins = options.get_int(option_section, "n_theta_bins", 20)
                theta_min = options.get_double(option_section, "theta_min", 1.) 
                theta_max = options.get_double(option_section, "theta_max", 300.)
                theta_min = arcmin_to_radians(theta_min)
                theta_max = arcmin_to_radians(theta_max)
                theta_edges = np.logspace(np.log10(theta_min), np.log10(theta_max), n_theta_bins + 1)
        print('theta_edges=',theta_edges)
        if theta_from_dat:
            theta = readtheta(filename, xi_type_2pt, theta_type = 'centers', desired_units = 'rad')
        else:
            if two_thirds_midpoint:
                theta = (2./3.) * (theta_edges[1:]**3 - theta_edges[:-1]**3) / (theta_edges[1:]**2 - theta_edges[:-1]**2)
                # this is currently the default
            else:
                theta = (theta_edges[:-1]*theta_edges[1:])**0.5
                # geometric mean
            # these are the values passed to the block and are needed for plotting or,
            # in some cases, for determining scale cuts.
    else:
        print("*** Using Legendre coefficients at individual theta values ***")
        if theta_from_dat:
            print("*** Reading in central theta values from datafile ***")
            theta = readtheta(filename, xi_type_2pt, theta_type = 'centers', desired_units = 'arcmin')
            theta = arcmin_to_radians(theta)
        elif options.has_value(option_section, "theta"):
            print('*** Note: Using `theta` values specified in ini file ***')
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

    # setup precompute functions and I/O sections

    if xi_type in ["22", "22+", "22-"]:
        if bin_avg:
            precomp_func = get_legfactors_22_binav
        else:    
            precomp_func = get_legfactors_22
        if not cl_section:
            cl_section = "shear_cl"
        if output_section == "":
            output_section = (cl_section.replace("cl","xi_plus"), cl_section.replace("cl","xi_minus"))
        else:
            output_section = (output_section + "_plus", output_section + "_minus")
    elif xi_type == "00":
        if bin_avg:
            precomp_func = get_legfactors_00_binav
        else:    
            precomp_func = get_legfactors_00
        if not cl_section:
            cl_section = "galaxy_cl"
        if output_section == "":
            output_section = cl_section.replace("cl", "xi")
    elif xi_type in ["02", "02+"]:
        if bin_avg:
            precomp_func = get_legfactors_02_binav
        else:    
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
    if bin_avg:
        legfacs = precomp_func(np.arange(ell_max + 1), theta_edges)
    else:
        legfacs = precomp_func(np.arange(ell_max + 1), theta)
    if high_l_filter>0:
        if isinstance(legfacs, tuple):
            legfacs = ( apply_filter( ell_max, high_l_filter, legfacs[0] ), 
            apply_filter( ell_max, high_l_filter, legfacs[1] ) )
        else:
            legfacs = apply_filter( ell_max, high_l_filter, legfacs )
    return theta, theta_edges, ell_max, legfacs, cl_section, output_section, save_name, bin_avg

def execute(block, config):

    thetas, theta_edges, ell_max, legfacs, cl_section, output_section, save_name, bin_avg = config
    print('thetas in execute are ',thetas*180*60/np.pi)

    ell = block[cl_section, "ell"]

    nbina, nbinb = block[cl_section, 'nbin_a'], block[cl_section, 'nbin_b']

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
        block[o, "theta_edges"] = theta_edges
        block[o, "bin_avg"] = bin_avg
        block.put_metadata(o, "theta", "unit", "radians")
    return 0

def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness.  The joy of python.
    return 0
