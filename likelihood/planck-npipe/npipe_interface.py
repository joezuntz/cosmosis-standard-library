from cosmosis.datablock import option_section
try:
    from planckpr4lensing.iswlens_jtliks.lik import cobaya_jtlik
except ImportError:
    raise ImportError("Please install the planckpr4lensing with: pip install git+https://github.com/carronj/planck_PR4_lensing")
import numpy as np


def setup(options):
    mode = options[option_section, "mode"]
    mode_choices = [
        'pp',
        'pt',
        'pt_pp',
        'tt',
        'tt_pp',
        'tt_pp_mtt',
        'tt_pt',
        'tt_pt_mtt',
        'tt_pt_pp',
        'tt_pt_pp_mpp',
        'tt_pt_pp_mtt',
        'tt_pt_pp_mtt_mpp']
    s = ", ".join(mode_choices)
    if mode not in mode_choices:
        raise ValueError(f"Mode parameter in npipe must be one of: {s}, not {mode}")
    calculator = cobaya_jtlik()
    calculator.initialize()
    chi2_method = getattr(calculator, f'get_chi2_{mode}')
    return calculator, chi2_method, mode



def execute(block, config):
    calculator, chi2_method, mode = config

    c_ells = {}
    ell = block['cmb_cl', 'ell']
    A = block['planck', 'a_planck']

    # Not sure if pt currently calculated in camb module - need to check. might be.
    # be careful with normalizations.
    c_ells['tt'] = block['cmb_cl', 'tt']
    if 'pp' in mode:
        c_ells['pp'] = block['cmb_cl', 'pp']
    if 'pt' in mode:
        c_ells['pt'] = block['cmb_cl', 'pt']
    if 'tt_pt' in mode:
        c_ells['ee'] = block['cmb_cl', 'ee']
        c_ells['te'] = block['cmb_cl', 'te']

    # npipe wants to start at zero.
    if ell[0] == 2:
        ell = np.concatenate([[0, 1], ell])
        for name, c_ell in c_ells.items():
            c_ells[name] = np.concatenate([[0.0, 0.0], c_ell])

    # Convert from D_ell to C_ell
    factor = ell * (ell + 1) / 2 / np.pi
    # avoid divide-by-zero
    factor[:2] = 1

    # Scale by calibration parameter as in this line:
    factor *= A**2

    for cl in c_ells.values():
        cl /= factor

    # Compute and save likelihood
    like = chi2_method(c_ells)    
    block['likelihoods', 'npipe_like'] = -like / 2

    return 0
