"""
Author: Sunao Sugiyama
Last edit: 2023.05.16

This does almost the same procedure as structure/minimalbias/minimalbias_interface.py: 
here we apply measurement corrections to the magnification bias term of dSigma.
"""
from cosmosis.datablock import names, option_section
from magnification_phys import magnification_class
import numpy as np
import os, sys
dirname = os.path.split(__file__)[0]
twopoint_path = os.path.join(dirname,"..","..","likelihood","2pt")
sys.path.insert(0, twopoint_path)
import twopoint
import time

def read_rp(filename, xi_type_2pt, theta_type = 'centers', desired_units = 'Mpc/h'):
    """
    Read the projected radial bins from data file.
    """
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
        raise ValueError("Unknown theta_type in read_theta")
        
        
def setup(options):
    """
    Defines the setups to calcurate the magnification bias term, and 
    instantiate necessary quantities.
    
    Now this module adds the magnification bias term only to g-g lensing signal, 
    but not to the galaxy clustering wp. (This is because that contribution is 
    small especially given the covariance of HSC weak lensing analysis.) 
    However, I (SS) leave some lines for wp for the future purpose.
    """
    
    # lens redshift: Note that dSigma and wp are independent of source redshift.
    redshifts = options.get_double_array_1d(option_section, 'redshifts')
    
    # Choice whether we average the signal within finite radial bin.
    bin_avg = options.get_bool(option_section, "bin_avg", False)
    
    # radial bins in reference cosmology for dSigma and wp
    theta_file = options.get_string(option_section, "theta_file", "")
    rp_ds = read_rp(theta_file, 'dsigma', theta_type = 'edges' if bin_avg else 'centers')
    rp_wp = read_rp(theta_file, 'wp'    , theta_type = 'edges' if bin_avg else 'centers')
    
    # Logrithmic bin width
    dlnrp_ds = np.log(rp_ds[1]/rp_ds[0])
    dlnrp_wp = np.log(rp_wp[1]/rp_wp[0])
    
    # Name of samples used for lens and source: this shoud be read by one of the number_density modules
    name = options.get_string(option_section, "lens_sample", "")
    section_lens = "NZ_" + name.upper()
    name = options.get_string(option_section, "source_sample", "")
    section_srce = "NZ_" + name.upper()
    
    # initialize magnification bias module
    config = {'verbose': options.get_bool(option_section, "verbose", True)}
    mag = magnification_class(config)
    
    # projection length
    pimax = options.get_double(option_section, 'pimax', default=100.0)
    
    # output section name
    output_section_wp = options.get_string(option_section, "output_section_name_wp", "")
    output_section_ds = options.get_string(option_section, "output_section_name_ds", "")
    
    if output_section_wp == "":
        output_section_wp = 'galaxy_wp'
    if output_section_ds == "":
        output_section_ds = 'galaxy_shear_dsigma'
    
    return redshifts, rp_ds, rp_wp, dlnrp_ds, dlnrp_wp, mag, section_lens, section_srce, pimax, output_section_wp, output_section_ds
    
    
def execute(block, config):
    """
    Name the section of magnification biass parameter as "magnification" 
    so that this modules can read the parameter.
    """
    time_start = time.time()
    
    # Get configuration
    redshifts, rp_ds, rp_wp, dlnrp_ds, dlnrp_wp, mag, section_lens, section_srce, pimax, output_section_wp, output_section_ds = config
    
    # number of redshift bins for lens and source.
    nbin_lens = block[section_lens, 'nbin']
    nbin_srce = block[section_srce, 'nbin']
    
    # power spectra
    z_nlin, k_nlin, P_nlin = block.get_grid(names.matter_power_nl, "z", "k_h", "P_k")
    mag.set_pk_nlin_data(z_nlin, k_nlin, P_nlin)
    
    # comoving distance
    h0  = block[names.cosmological_parameters, 'h0']
    z   = block[names.distances, "z"]
    chi = block[names.distances, "D_C"] * h0
    mag.set_z2chi(z, chi)
    
    # iterate over all lens-source bin pair
    for i in range(nbin_lens):
        # Get representative redshift of i-th lens bin
        zl = redshifts[i]
        
        # Linear galayx bias
        alpha  = block['magnification_bias_parameters', 'alpha_{0}'.format(i+1)]
        Omm = block[names.cosmological_parameters, "Omega_m"]
        mag.set_param(alpha, Omm)
        
        # lens redshift distribution
        mag.set_nz_lens(block[section_lens, "z"], block[section_lens, "bin_{0}".format(i+1)])
        
        # correction factors 
        f_rp = block.get_double('meascorr_rp', 'bin_{0}'.format(i+1), 1.0)
        f_wp = block.get_double('meascorr_wp', 'bin_{0}'.format(i+1), 1.0) # multiplied to pimax
        
        # update wp
        block[output_section_wp+'_gg', 'bin_{0}'.format(i+1)] = block[output_section_wp, 'bin_{0}'.format(i+1)]
        block[output_section_wp, 'bin_{0}'.format(i+1)] += 0.0 #mag.get_wp_mag(z, rp ,dlnrp,f_wp*pimax)
        
        # apply correction to dSigma overall amplitude
        for j in range(nbin_srce):
            mag.set_nz_source(block[section_srce, "z"], block[section_srce, "bin_{0}".format(j+1)])
            f_ds = block.get_double('meascorr_ds', 'bin_{0}_{1}'.format(i+1,j+1), 1.0)
            block[output_section_ds+'_gG', 'bin_{0}_{1}'.format(i+1,j+1)] = block[output_section_ds, 'bin_{0}_{1}'.format(i+1,j+1)]
            block[output_section_ds, 'bin_{0}_{1}'.format(i+1,j+1)] += f_ds * mag.get_ds_mag(zl, f_rp*rp_ds, dlnrp_ds)
    
    # record elapsed time.
    block['elapsed_time', 'magnification'] = time.time() - time_start
    block.put_metadata('elapsed_time', "magnification", "unit", "seconds")
    
    return 0