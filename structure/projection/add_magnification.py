"""
add the magnification terms to the shear-gal and gal-gal Cls
based on the add_intrinsic module
"""
from __future__ import print_function
from builtins import range
from cosmosis.datablock import option_section, names
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

def setup(options):
    do_galaxy_galaxy = options.get_bool(option_section, "galaxy-galaxy", True)
    do_galaxy_shear = options.get_bool(option_section, "galaxy-shear", True)
    include_intrinsic = options.get_bool(option_section, "include_intrinsic", False)
    print()
    print("The add_magnification module will try to combine magnification terms into")
    if do_galaxy_galaxy and do_galaxy_shear:
        print("both the galaxy-galaxy and galaxy-shear spectra.")
    elif do_galaxy_galaxy:
        print("only the galaxy-galaxy spectra.")
    elif do_galaxy_shear:
        print("only the galaxy-shear.")
    else:
        print("... actually not into anything. You set galaxy-galaxy=F and galaxy-shear=F")
        print("Ths module will not do anything in this configuration")
    print()
    return do_galaxy_galaxy, do_galaxy_shear, include_intrinsic

def execute(block, config):
    do_galaxy_galaxy, do_galaxy_shear, include_intrinsic = config
    if do_galaxy_galaxy:
        nbin_pos = block[names.galaxy_cl, 'nbin']
    elif do_galaxy_shear:
        nbin_pos = block["galaxy_shear_cl", 'nbin_a']
    if do_galaxy_shear:
        nbin_shear = block["galaxy_shear_cl", 'nbin_b']
    if do_galaxy_galaxy:
        # for galaxy_galaxy, we're replacing 'galaxy_cl' (the gg term) with gg+gm+mg+mm
        # so in case useful, save the gg term to galaxy_cl_gg.
        ells = block[names.galaxy_cl, 'ell']
        block["galaxy_cl_gg", 'ell'] = ells
        gal_mag_ells = block["galaxy_magnification_cl", 'ell']
        mag_mag_ells = block["magnification_cl", "ell"]
        resample_gal_mag = False
        resample_mag_mag = False
        if len(gal_mag_ells) != len(ells):
            resample_gal_mag = True
        elif not np.allclose(gal_mag_ells, ells):
            resample_gal_mag = True
        if len(mag_mag_ells) != len(ells):
            resample_mag_mag = True
        elif not np.allclose(mag_mag_ells, ells):
            resample_mag_mag = True
        #Get auto_only - if True, only auto-correlations for galaxy_cl
        #were computed.
        auto_only = block.get_bool(names.galaxy_cl, "auto_only", False)
        for i in range(nbin_pos):
            for j in range(i + 1):
                if (i!=j and auto_only):
                    continue
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                bin_ji = 'bin_{1}_{0}'.format(i + 1, j + 1)
                block["galaxy_cl_gg", bin_ij] = block[names.galaxy_cl, bin_ij]
                gal_mag_ij_orig = block["galaxy_magnification_cl", bin_ij]
                gal_mag_ji_orig = block["galaxy_magnification_cl", bin_ji]
                mag_mag_orig = block["magnification_cl", bin_ij]
                if resample_gal_mag:
                    gal_mag_ij = np.zeros_like(ells)
                    gal_mag_ji = np.zeros_like(ells)
                    mag_mag = np.zeros_like(ells)
                    gal_mag_ij[ells>0] = InterpolatedUnivariateSpline(np.log(gal_mag_ells), 
                        gal_mag_ij_orig)(np.log(ells[ells>0]))
                    gal_mag_ji[ells>0] = InterpolatedUnivariateSpline(np.log(gal_mag_ells), 
                        gal_mag_ji_orig)(np.log(ells[ells>0]))  
                else:
                    gal_mag_ij = gal_mag_ij_orig
                    gal_mag_ji = gal_mag_ji_orig
                if resample_mag_mag:
                    mag_mag[ells>0] = InterpolatedUnivariateSpline(np.log(gal_mag_ells), 
                        gal_mag_ji_orig)(np.log(ells[ells>0]))   
                else:
                    mag_mag = mag_mag_orig    
                block[names.galaxy_cl, bin_ij] += (
                    gal_mag_ij + gal_mag_ji
                    + mag_mag   #mm
                )
    #if include_intrinsic is True, we're replacing gG+gI with gG+gI+mG+mI       
    #else,                         we're replacing gG    with gG+mG 
    if do_galaxy_shear:
        for i in range(nbin_pos):
            for j in range(nbin_shear):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                block["galaxy_shear_cl_gG", bin_ij] = block["galaxy_shear_cl", bin_ij]
                block["galaxy_shear_cl", bin_ij] += (
                    block["magnification_shear_cl", bin_ij]
                    )
                if include_intrinsic == True:
                    block["galaxy_shear_cl", bin_ij] += (
                        block["magnification_intrinsic_cl", bin_ij]
                        )
    return 0