from cosmosis.datablock import option_section, names
from scipy.interpolate import interp1d
from numpy import zeros_like, zeros
import numpy as np

# this module multiplies the "target_section" by the nonlinear boost (matter_power_nl/matter_power_lin). It assumes the target nonlinearity is traced by the matter perturbation field nonlinearity.

def setup(options):
    lin = options.get_string(option_section, 'linear_section', names.matter_power_lin)
    nonlin = options.get_string(option_section, 'nonlinear_section', names.matter_power_nl)
    target = options.get_string(option_section, 'target_section', 'weyl_curvature_spectrum')
    return lin, nonlin, target

def execute(block, config):
    lin, nonlin, target = config

    # Need to get things at the same wavenumbers for ratio
    z_lin,k_lin,p_k_lin = block.get_grid(lin, "z","k_h", "P_k")
    z_nl,k_nl,p_k_nl = block.get_grid(nonlin, "z", "k_h","P_k")
    z_t,k_t,p_k_t = block.get_grid(target, "z","k_h", "P_k")
    
    if not (np.allclose(z_nl, z_lin) & np.allclose(z_t, z_lin)):
        raise ValueError("NLFactor module requires the same z values for lin, nl matter power and the target power.")

    p_k_lin_at_nlk = zeros_like(p_k_nl)
    # assuming same z array

    for zi in range(len(z_nl)):
        p_k_lin_at_nlk[zi,:] = interp1d(k_lin, p_k_lin[zi,:])(k_nl)

    R = p_k_nl / p_k_lin_at_nlk

    # The nonlin pk doesn't go to as small k, so we have to truncate the target k_h for the nl case.
    if min(k_t)<min(k_nl):
        mink = next(j[0] for j in enumerate(k_t) if j[1]>min(k_nl))
    else:
        mink = 0
    if max(k_t)>max(k_nl):
        maxk = next(j[0] for j in enumerate(k_t) if j[1]>max(k_nl))
    else:
        maxk=len(k_t)
    ratio_at_wk = zeros((len(z_nl), len(k_t[mink:maxk])))
    for zi in range(0, len(z_nl)):
        ratio_at_wk[zi,:] = interp1d(k_nl, R[zi,:])(k_t[mink:maxk])
    
    # Apply ratio to target spectrum
    weyl_nl = zeros_like(ratio_at_wk)
    for zi in range(0, len(z_nl)):
        weyl_nl[zi,:] = ratio_at_wk[zi,:]*p_k_t[zi,:][mink:maxk]   
    
    block.put_grid(target+"_nl","z",z_nl,"k_h",k_t[mink:maxk],"P_k",weyl_nl)

    return 0
