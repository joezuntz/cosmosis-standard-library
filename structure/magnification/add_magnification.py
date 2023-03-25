"""
add the magnification terms to the shear-gal and gal-gal Cls
based on the add_intrinsic module
"""
from cosmosis.datablock import option_section, names
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

def setup(options):
    do_galaxy_galaxy = options.get_bool(option_section, "galaxy-galaxy", True)
    do_galaxy_shear = options.get_bool(option_section, "galaxy-shear", True)
    do_galaxy_cmbkappa = options.get_bool(option_section, "galaxy-cmbkappa", False)
    include_intrinsic = options.get_bool(option_section, "include_intrinsic", False)
    print()
    print("The add_magnification module will try to combine magnification terms into the following spectra:")
    if do_galaxy_shear:
        print("galaxy-shear")        
    if do_galaxy_galaxy:
        print("galaxy-galaxy")
    if do_galaxy_cmbkappa:
        print("galaxy-cmbkappa")
    if (not do_galaxy_shear*do_galaxy_galaxy*do_galaxy_cmbkappa):
        print("... actually not into anything. You set galaxy-galaxy=F and galaxy-shear=F (and possibly galaxy-cmbkappa = False)")
        print("Ths module will not do anything in this configuration")
    print()
    return do_galaxy_galaxy, do_galaxy_shear, do_galaxy_cmbkappa, include_intrinsic

def execute(block, config):
    do_galaxy_galaxy, do_galaxy_shear, do_galaxy_cmbkappa, include_intrinsic = config

    # Get the number of bins for the different spectra.
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

        # The different C_ell might conceivably be defined on different grids.
        # If so we will need to interpolate
        gal_mag_ells = block["galaxy_magnification_cl", 'ell']
        mag_mag_ells = block["magnification_cl", "ell"]
        resample_gal_mag = (len(gal_mag_ells) != len(ells)) or not np.allclose(gal_mag_ells, ells)
        resample_mag_mag = (len(mag_mag_ells) != len(ells)) or not np.allclose(mag_mag_ells, ells)


        if resample_gal_mag or resample_mag_mag:
            index = ells > 0
            log_ell = np.log(ells[index])
            log_ell_gm = np.log(gal_mag_ells)
            log_ell_mm = np.log(mag_mag_ells)

        #Get auto_only - if True, only auto-correlations for galaxy_cl
        #were computed.
        auto_only = block.get_bool(names.galaxy_cl, "auto_only", False)
        for i in range(nbin_pos):
            for j in range(i + 1):
                if (i!=j and auto_only):
                    continue
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                bin_ji = 'bin_{1}_{0}'.format(i + 1, j + 1)

                # Load the C_ell for this bin pair, for the different spectra
                block["galaxy_cl_gg", bin_ij] = block[names.galaxy_cl, bin_ij]
                gal_mag_ij_orig = block["galaxy_magnification_cl", bin_ij]
                gal_mag_ji_orig = block["galaxy_magnification_cl", bin_ji]
                mag_mag_orig = block["magnification_cl", bin_ij]


                # Optionally resample the GM term, if needed
                if resample_gal_mag:
                    gal_mag_ij = np.zeros_like(ells)
                    gal_mag_ji = np.zeros_like(ells)
                    mag_mag = np.zeros_like(ells)
                    gal_mag_ij[index] = InterpolatedUnivariateSpline(log_ell_gm,
                        gal_mag_ij_orig)(log_ell)
                    gal_mag_ji[index] = InterpolatedUnivariateSpline(log_ell_gm,
                        gal_mag_ji_orig)(log_ell)
                else:
                    gal_mag_ij = gal_mag_ij_orig
                    gal_mag_ji = gal_mag_ji_orig

                # Optionally resample the MM term, if needed
                if resample_mag_mag:
                    mag_mag[index] = InterpolatedUnivariateSpline(log_ell_mm,
                        mag_mag_orig)(log_ell)
                else:
                    mag_mag = mag_mag_orig    

                # Finally, combine all the terms together
                block[names.galaxy_cl, bin_ij] += gal_mag_ij + gal_mag_ji + mag_mag

    #if include_intrinsic is True, we're replacing gG+gI with gG+gI+mG+mI       
    #else,                         we're replacing gG    with gG+mG 
    if do_galaxy_shear:
        block["galaxy_shear_cl_gG", "ell"] = block["galaxy_shear_cl", "ell"]
        for i in range(nbin_pos):
            for j in range(nbin_shear):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                block["galaxy_shear_cl_gG", bin_ij] = block["galaxy_shear_cl", bin_ij]
                block["galaxy_shear_cl", bin_ij] += (
                    block["magnification_shear_cl", bin_ij]
                    )
                if include_intrinsic:
                    block["galaxy_shear_cl", bin_ij] += (
                        block["magnification_intrinsic_cl", bin_ij]
                        )

    if do_galaxy_cmbkappa:
        for i in range(nbin_pos):
            bin_i = 'bin_{0}_1'.format(i + 1)
            block["galaxy_cmbkappa_cl_gK", bin_i] = block["galaxy_cmbkappa_cl", bin_i]
            block["galaxy_cmbkappa_cl", bin_i] += (
                block["magnification_cmbkappa_cl", bin_i]
                )
                    
    return 0
