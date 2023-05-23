"""
Author: Sunao Sugiyama
Last edit: 2023.05.16

This module computes the 2x2pt observables in physical scale. This choice is used in HSC weak lensing analysis.
The observables are the galaxy-galaxy lensign, dSigma(r_p), and the projected correlation function, w_p(r_p).
This module returns the observables after measureemnt corrections:
    1. correction on projected radial bin length r_p -> f * r_p, where f is the ratio 
        of comoving distance at the lens redshift and correct for the difference of cosmologies 
        in measurement and modeling.
    2. correction on projection length of w_p, pimax -> f * pimax, where f is the ratio of 
        Hubble parameter at the lens redshift and correct for the difference of cosmologies 
        in measurement and modeling.
    3. correction on overall dSigma amplitude, dSigma -> f * dSigma, where f is the ratio of 
        Sigma_crit values averaged over all lens-source pairs, and corrects for the differene of 
        cosmologies and photometric redshift in measurement and modeling.
This module also take account of the averaging effect within the finite radial bin.

The radial bins are read in this module at the time of setup, and the correction factors introduced above are all computed in a module of 'shear/meascorr/meascorr_interface.py', and the results are saved in data block.
"""
from cosmosis.datablock import names, option_section
from minimalbias import minimalbias_class
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
    
    # initialize minimal bias module
    mb_config  = {'do_Kaiser': options.get_bool(option_section, 'do_Kaiser', default=True), 
                  'verbose'  : options.get_bool(option_section, 'verbose', default=True)}
    mb = minimalbias_class(mb_config)
    
    # projection length
    pimax = options.get_double(option_section, 'pimax', default=100.0)
    
    # output section name
    output_section_wp = options.get_string(option_section, "output_section_name_wp", "")
    output_section_ds = options.get_string(option_section, "output_section_name_ds", "")
    
    if output_section_wp == "":
        output_section_wp = 'galaxy_wp'
    if output_section_ds == "":
        output_section_ds = 'galaxy_shear_dsigma'
    
    return redshifts, rp_ds, rp_wp, dlnrp_ds, dlnrp_wp, mb, section_lens, section_srce, pimax, output_section_wp, output_section_ds
    
    
def execute(block, config):
    time_start = time.time()
    
    # Get configuration
    redshifts, rp_ds, rp_wp, dlnrp_ds, dlnrp_wp, mb, section_lens, section_srce, pimax, output_section_wp, output_section_ds = config
    
    # number of redshift bins for lens and source.
    nbin_lens = block[section_lens, 'nbin']
    nbin_srce = block[section_srce, 'nbin']
    
    # power spectra
    z_lin , k_lin , P_lin  = block.get_grid(names.matter_power_lin, "z", "k_h", "P_k")
    z_nlin, k_nlin, P_nlin = block.get_grid(names.matter_power_nl, "z", "k_h", "P_k")
    mb.set_pk_lin_data(z_lin , k_lin , P_lin)
    mb.set_pk_nlin_data(z_nlin, k_nlin, P_nlin)
    
    # z & chi relation
    z   = block[names.growth_parameters, "z"]
    f   = block[names.growth_parameters, "f_z"]
    mb.set_z2f(z, f)
    
    # iterate over all lens-source bin pair
    for i in range(nbin_lens):
        # Get representative redshift of i-th lens bin
        zl = redshifts[i]
        
        # Linear galayx bias
        b1  = block['galaxy_bias_parameters', 'b1_{0}'.format(i+1)]
        Omm = block[names.cosmological_parameters, "Omega_m"]
        mb.set_param(Omm, b1)
        
        # correction factors 
        f_rp = block.get_double('meascorr_rp', 'bin_{0}'.format(i+1), 1.0)
        f_wp = block.get_double('meascorr_wp', 'bin_{0}'.format(i+1), 1.0) # multiplied to pimax
        
        # update wp
        block[output_section_wp, 'bin_{0}'.format(i+1)] = mb.get_wp(zl, f_rp*rp_wp, dlnrp_wp, f_wp*pimax)
        
        # apply correction to dSigma overall amplitude
        for j in range(nbin_srce):
            f_ds = block.get_double('meascorr_ds', 'bin_{0}_{1}'.format(i+1,j+1), 1.0)
            block[output_section_ds, 'bin_{0}_{1}'.format(i+1,j+1)] = f_ds * mb.get_ds(zl, f_rp*rp_ds, dlnrp_ds)
    
    # save radial bin in reference cosmology (i.e. without measurement correction)
    block[output_section_wp, 'rp'] = rp_wp
    block.put_metadata(output_section_wp, "rp", "unit", "Mpc/h")
    block[output_section_ds, 'rp'] = rp_ds
    block.put_metadata(output_section_ds, "rp", "unit", "Mpc/h")
    
    # record elapsed time.
    block['elapsed_time', 'minimalbias'] = time.time() - time_start
    block.put_metadata('elapsed_time', "minimalbias", "unit", "seconds")
    
    return 0