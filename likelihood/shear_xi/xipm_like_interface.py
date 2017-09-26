from __future__ import print_function
# Module to compute a general xi+/- likelihood.
# Signal should be in ascii file with columns: theta xip_11, xip_12,...,xip_22,xip_23,....xim_11 etc.
# Covariance in .npy or ascii file with same ordering.

from builtins import range
from cosmosis.datablock import option_section, names as section_names
import xipm_like
import numpy as np
import os


def setup(options):
    sec = option_section
    covmat_file = options.get_string(sec, 'covariance_file')
    data_file = options.get_string(sec, 'data_file')
    n_z_bins = options.get_int(sec, 'n_z_bins')
    cov_num_rlzn = options.get_int(sec, 'cov_num_rlzn', default=0)
    pmc = options.get_bool(sec, 'pmc', default=False)
    if cov_num_rlzn == 0:
        cov_num_rlzn = None
    # This allows you do load a theory signal from cosmosis output
    # for debugging/forecasting
    theory_signal_dir = options.get_string(
        sec, 'theory_signal_dir', default='')
    if theory_signal_dir == '':
        theory_signal_dir = None
    if theory_signal_dir is not None:
        print('using theory signal from %s' % theory_signal_dir)
        import scipy.interpolate as interp

    plus_only = options.get_bool(sec, 'plus_only', default=False)
    theta_mins = options.get_string(sec, 'theta_mins', default='')
    if theta_mins == '':
        theta_mins = None
    else:
        theta_mins = eval(theta_mins)

    print('Minimum theta values = ', theta_mins)
    theta_maxs = options.get_string(sec, 'theta_maxs', default='')
    if theta_maxs == '':
        theta_maxs = None
    else:
        theta_maxs = eval(theta_maxs)
    print('Maximum theta values = ', theta_maxs)

    # create likelihood calculator
    # loads named files and prepares itself
    """
    try:
        covmat=np.load(covmat_file)
    except IOError:
        try:
            covmat=np.loadtxt(covmat_file)
        except IOError:
            print "couldn't load covariance matrix as pickle or ascii file"
            return 1
    """
    data = np.loadtxt(data_file).T
    theta_data = data[0]
    xipm_data = {}
    xipm_data['theta'] = theta_data

    k = 0
    n_z_bin_pairs = n_z_bins * (n_z_bins + 1) / 2
    for i in range(1, n_z_bins + 1):
        for j in range(i, n_z_bins + 1):
            bin_comb = '%d_%d' % (j, i)
            if theory_signal_dir == None:
                xipm_data['xip_' + bin_comb] = data[k + 1]
                xipm_data['xim_' + bin_comb] = data[n_z_bin_pairs + k + 1]
            else:
                theory_theta = 60. * \
                    np.degrees(np.loadtxt(os.path.join(
                        theory_signal_dir, 'theta.txt')).T)
                theory_xip = np.loadtxt(os.path.join(
                    theory_signal_dir, 'xiplus_%d_%d.txt' % (j, i))).T
                theory_xim = np.loadtxt(os.path.join(
                    theory_signal_dir, 'ximinus_%d_%d.txt' % (j, i))).T
                xip = np.exp(interp.interp1d(np.log(theory_theta),
                                             np.log(theory_xip))(np.log(theta_data)))
                xim = np.exp(interp.interp1d(np.log(theory_theta),
                                             np.log(theory_xim))(np.log(theta_data)))
                xipm_data['xip_' + bin_comb] = xip
                xipm_data['xim_' + bin_comb] = xim
            k += 1

    calculator = xipm_like.XipmLikelihood(covmat_file,
                                          xipm_data, n_z_bins, cov_num_rlzn=cov_num_rlzn, theta_mins=theta_mins, theta_maxs=theta_maxs, plus_only=plus_only, pmc=pmc)

    # pass back config to the
    return calculator, n_z_bins


def execute(block, config):
    calculator, n_z_bins = config

    # Get theta for the sample values.
    # We need theta twice because our CFHTLens
    # code wants xminus and xplus
    section = section_names.shear_xi
    theta = block[section, "theta"]
    theta = np.concatenate((theta, theta))

    # Get the xi(theta) for these samples, for each pair of bins.
    # The likelihood calculator wants a big dictionary
    xi_theory = {}
    for i in range(1, n_z_bins + 1):
        for j in range(i, n_z_bins + 1):
            name = 'xiplus_%d_%d' % (j, i)
            xiplus = block[section, name]
            name = 'ximinus_%d_%d' % (j, i)
            ximinus = block[section, name]
            xi = np.concatenate((xiplus, ximinus))
            xi_theory[(i, j)] = (theta, xi)

    # Calculate the likelihood

    like, theory_vector, data_vector, inv_cov, cov, sim = calculator(xi_theory)

    # save the result
    section = section_names.likelihoods
    block[section, "xipm_like"] = like

    # Also save the data vector
    block[section_names.data_vector, "xipm_theory"] = theory_vector
    block[section_names.data_vector, "xipm_data"] = data_vector
    block[section_names.data_vector, "xipm_inverse_covariance"] = inv_cov
    block[section_names.data_vector, "xipm_covariance"] = cov
    block[section_names.data_vector, "xipm_simulation"] = sim

    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness.  The joy of python.
    return 0
