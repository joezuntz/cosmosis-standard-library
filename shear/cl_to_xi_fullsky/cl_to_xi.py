from builtins import object
import numpy as np
import scipy.interpolate as interp
from legendre import *
from numpy import pi, cos, sin, sqrt


class SpectrumInterp(object):
    """class for 1d interpolation of spectra"""

    def __init__(self, angle, spec, kind='cubic'):
        if np.all(spec > 0):
            self.interp_func = interp.interp1d(np.log(angle), np.log(spec), bounds_error=False,
                                               kind=kind)
            self.interp_type = 'loglog'
        elif np.all(spec < 0):
            self.interp_func = interp.interp1d(np.log(angle), np.log(-spec), bounds_error=False,
                                               kind=kind)
            self.interp_type = 'minus_loglog'
        else:
            self.interp_func = interp.interp1d(np.log(angle), spec, bounds_error=False, fill_value=0.,
                                               kind=kind)
            self.interp_type = "log_ang"

    def __call__(self, angle):
        if self.interp_type == 'loglog':
            spec = np.exp(self.interp_func(np.log(angle)))
        elif self.interp_type == 'minus_loglog':
            spec = -np.exp(self.interp_func(np.log(angle)))
        else:
            assert self.interp_type == "log_ang"
            spec = self.interp_func(np.log(angle))
        spec[np.isnan(spec)] = 0.
        return spec

def radians_to_arcmin(r):
    return np.degrees(r) * 60

def arcmin_to_radians(a):
    return np.radians(a / 60.)

def cl_to_xi_to_block(block, output_section, name,
                          cl_interp, thetas, legfacs):
    if not isinstance(output_section, tuple):
        output_section = (output_section,)
        legfacs = (legfacs,)
    for (o,leg) in zip( output_section, legfacs ):
        xis = np.zeros_like(thetas)
        ell_max = leg.shape[1] - 1
        try:
            assert len(cl_interp) == leg.shape[1]
        except TypeError:
            ells = np.arange(ell_max + 1)
            cl_interp = cl_interp(ells)

        for it, t in enumerate(thetas):
            xis[it] = np.sum(cl_interp * leg[it])
        block[o, name] = xis

def cl_to_xi_precomp_00_02(cl_input, thetas, legfacs):
    """This function works for spin (0,0) e.g. w(theta)
    and spin 0,2 i.e gamma_t."""
    xis = np.zeros_like(thetas)
    ell_max = legfacs.shape[1] - 1
    try:
        assert len(cl_input) == legfacs.shape[1]
    except TypeError:
        ells = np.arange(ell_max + 1)
        cl_input = cl_input(ells)

    for it, t in enumerate(thetas):
        xis[it] = np.sum(cl_input * legfacs[it])
    return xis

def cl_to_xi_plus_and_minus_precomp(cl_input, thetas, G_plus_minus_pre):
    """Compute xi+/- given:
    cl_input: Interplation function for cl
    ell_max: maximum ell (will start from ell=2)
    thetas: real space angles to calculate xi+/-(theta)
    G_plus_minus_pre: tuple containing G+ and G- from precomp_GpGm. See eqns. 2.26,2.27 
    from astro-ph/9611125v1"""
    G_plus_pre, G_minus_pre = G_plus_minus_pre
    ell_max = G_plus_pre.shape[1] - 1
    ells = np.arange(ell_max + 1)
    N_ell = get_N_ell(ells)
    N_ell[:2] = 0.
    assert G_plus_pre.shape[0] == G_minus_pre.shape[0] == len(thetas)
    assert G_plus_pre.shape[1] == G_minus_pre.shape[1] == ell_max + 1
    C_plus, C_cross = np.zeros_like(thetas), np.zeros_like(thetas)
    Cls = cl_input(ells)
    for i, theta in enumerate(thetas):
        # print theta, G_plus, G_minus
        C_plus[i] = np.sum(((2 * ells + 1) / 2. / pi) *
                           N_ell**2 * Cls * G_plus_pre[i]) / 2
        C_cross[i] = np.sum(((2 * ells + 1) / 2. / pi) *
                            N_ell**2 * Cls * G_minus_pre[i]) / 2

    xi_plus = C_plus + C_cross
    xi_minus = C_plus - C_cross

    return xi_plus, xi_minus


def save_xi_00_02(block, section, bin_i, bin_j, cl_input, thetas, legfactors):
    name = "bin_%d_%d" % (bin_i, bin_j)
    w = cl_to_xi_precomp_00_02(cl_input, thetas, legfactors)
    block[section, name] = w
    return w

def save_xi_22(block, section, bin_i, bin_j, cl_input, thetas, G_plus_minus_pre):
    xi_plus, xi_minus = cl_to_xi_plus_and_minus_precomp(
        cl_input, thetas, G_plus_minus_pre)
    block[section[0], "bin_%d_%d" % (bin_i, bin_j)] = xi_plus
    block[section[1], "bin_%d_%d" % (bin_i, bin_j)] = xi_minus
    return xi_plus, xi_minus
