from __future__ import print_function
from builtins import range
import sys
import os
import numpy as np
import cosmosis_py
from cosmosis.datablock import names, option_section
import cosmosis_py.section_names

sn_params = names.supernova_params
distances = names.distances
likes = names.likelihoods
options = option_section


run_dir = os.path.split(__file__)[0]  # detect path of this directory
# utils_dir=os.path.join(run_dir,"..","sn_utils") #To get path of utils dir
# sys.path.append(utils_dir) #Append utils dir to path

import kb_utils  # K. Barbary file reading utils

INT_SIG_DEFAULT = 0.1  # Currently not fitting for sig_int
int_sig = INT_SIG_DEFAULT


def setup(block):
    datadir = block.get_string(options, "dirname", default=run_dir)
    datafile = block.get_string(options, "filename", "SALT2.FITRES")
    datapath = os.path.join(datadir, datafile)

    meta, sntab = kb_utils.read_datafile(datapath,
                                         default_tablename='SN', output_ndarray=False)

    sntab = sntab['SN']

    def get(name):
        try:
            return np.array(sntab[name.lower()])
        except KeyError:
            return np.array(sntab[name])
    # Pull out the columns we need
    z_obs = get('Z')
    z_sig = get("ZERR")
    c_obs = get("c")
    sig_c = get('cERR')
    x1_obs = get('x1')
    sig_x1 = get('x1ERR')
    x0_obs = get('x0')
    sig_x0 = get('x0ERR')
    mb_obs = get("mB")
    sig_mb = get("mBERR")
    covx0x1 = get('COVx0x1')
    covx0c = get('COVx0c')
    covx1c = get('COVx1c')

    scalefac = -2.5 / (x0_obs * np.log(10.0))
    covmbx1 = covx0x1 * scalefac
    covmbc = covx0c * scalefac

    covmats = []
    ok_covmats = []
    nbad = 0

    for i_c in range(len(z_obs)):
        covmat = np.array([[sig_mb[i_c]**2, covmbx1[i_c], covmbc[i_c]],
                           [covmbx1[i_c], sig_x1[i_c]**2, covx1c[i_c]],
                           [covmbc[i_c], covx1c[i_c], sig_c[i_c]**2]])

        if (covmat.diagonal().min() > 0 and np.isfinite(covmat).all()
                and np.linalg.eigvals(covmat).min() > 0):
            ok_covmats.append(True)
            covmats.append(covmat)
        else:
            ok_covmats.append(False)
            nbad += 1
    if nbad > 0:
        print("Note: the SALT2 process can produce unphysical covariance matrices")
        print("which are not finite.  This happened for %d out of the %d supernovae in this sample." % (nbad, len(z_obs)))
        print("These will not be used.  (This is not an error, just an advisory note.)")

    ok_covmats = np.array(ok_covmats)
    z_obs = z_obs[ok_covmats]
    z_sig = z_sig[ok_covmats]
    c_obs = c_obs[ok_covmats]
    x1_obs = x1_obs[ok_covmats]
    mb_obs = mb_obs[ok_covmats]

    return z_obs, z_sig, c_obs, x1_obs, mb_obs, covmats


def likelihood(data_vec, z_model_table, mu_model_table, M0, alpha, beta):

    z_obs, z_sig, c_obs, x1_obs, mb_obs, covmats = data_vec

    # calculate mu_obs from data and M0, alpha, beta
    mu_obs = mb_obs - M0 + alpha * x1_obs - beta * c_obs
    # calculate true mu from look up table
    mu_theory = np.interp(z_obs, z_model_table, mu_model_table)
    # Make array of nuisance parameters
    psi = np.array([-1, alpha, beta])

    # Build up the vector of sigma&2 for each SN
    mu_sig_sq = np.zeros(len(z_obs))
    for i_c in range(len(z_obs)):
        fit_sig_sq = np.dot(np.dot(psi, covmats[i_c]), psi)
        mu_sig_sq[i_c] = fit_sig_sq + int_sig**2

    # Get overall error vector
    chisquare = (mu_obs - mu_theory)**2 / mu_sig_sq

    # Return log likelihood = -chi^2/2
    LogLike = -0.5 * chisquare.sum()
    return LogLike


def execute(block, config):
    # Get nuisance parameters
    M0 = block[sn_params, 'M0']
    alpha = block[sn_params, 'alpha']
    beta = block[sn_params, 'beta']

    # Get mu(z) for theory model
    z_model_table = block[distances, 'Z']
    mu_model_table = block[distances, 'MU']

    # Make sure z is increasing
    if z_model_table[1] < z_model_table[0]:
        z_model_table = z_model_table[::-1]
        mu_model_table = mu_model_table[::-1]

    # calculate the log-likelihood
    LogLike = likelihood(config, z_model_table,
                         mu_model_table, M0, alpha, beta)

    # Give a little warning about infinity and NaN errors
    if not np.isfinite(LogLike):
        sys.stderr.write("Non-finite LogLike in sn_SALT2chi_like\n")

    # Save the result
    block[likes, 'SN_LIKE'] = float(LogLike)

    # Signal success
    return 0
