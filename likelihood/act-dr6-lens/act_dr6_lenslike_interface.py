import act_dr6_lenslike
from cosmosis.datablock import option_section, names
import numpy as np

def setup(options):
    variant = options.get_string(option_section, 'variant', default='act_baseline')
    lens_only = options.get_bool(option_section, 'lens_only', default=False)
    like_corrections = options.get_bool(option_section, 'like_corrections', default=False)
    # lens_only use True if not combining with any primary CMB data
    # like_corrections should be False if lens_only is True

    # This dict will now have entries like `data_binned_clkk` (binned data vector), `cov`
    # (covariance matrix) and `binmat_act` (binning matrix to be applied to a theory
    # curve starting at ell=0).
    data_dict = act_dr6_lenslike.load_data(variant,lens_only=lens_only,like_corrections=like_corrections)
    return data_dict



def execute(block, config):
    data_dict = config

    # These are the CMB lensing convergence spectra (not potential or deflection)
    # as well as the TT, EE, TE, BB CMB spectra (needed for likelihood corrections)
    # in uK^2 units. All of these are C_ell (not D_ell), no ell or 2pi factors.
    ell = block[names.cmb_cl, 'ell']
    f1 = ell * (ell + 1) / 2 / np.pi
    cl_tt = block[names.cmb_cl, 'tt'] / f1
    cl_ee = block[names.cmb_cl, 'ee'] / f1
    cl_te = block[names.cmb_cl, 'te'] / f1
    cl_bb = block[names.cmb_cl, 'bb'] / f1

    # Slightly different normalization here
    f2 = ell * (ell + 1)
    cl_pp = block[names.cmb_cl, 'pp'] / f2

    cl_kk = act_dr6_lenslike.pp_to_kk(cl_pp, ell)

    # Then call
    lnlike = act_dr6_lenslike.generic_lnlike(data_dict,ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb)
    block[names.likelihoods, 'act_dr6_lens_like'] = lnlike

    return 0