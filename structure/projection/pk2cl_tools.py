#python code for P(k) to C(l) calculations including
#Limber and exact ("non-Limber"). If this is too slow
#should be straightforward to convert to c.
from __future__ import print_function
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.integrate import quad
import sys
from LOG_HT import fft_log
import time

inv_sqrt2pi = 1./np.sqrt(2*np.pi)

def limber_integral(kernel1, kernel2, pk_interp_logk, ells, chimin, chimax, dchi,
    method="trapz"):
    print("""Doing Limber integral with method %s between 
        chi_min: %.2e and chi_max: %.2e with step size %.2e"""%(method, chimin, chimax, dchi))
    assert chimin>0.
    c_ells, c_ell_errs = np.zeros_like(ells), np.zeros_like(ells)
    chi_vals = np.arange(chimin, chimax+dchi, dchi)
    kernel1_vals = kernel1(chi_vals)
    kernel2_vals = kernel2(chi_vals)
    k1k2 = kernel1_vals * kernel2_vals

    if method == "trapz":
        K_VALS = (ells[:, np.newaxis]+0.5) / chi_vals
        CHI_VALS = (chi_vals[:, np.newaxis]).T * np.ones((ells.shape[0], chi_vals.shape[0]))
        K1K2 = (k1k2[:, np.newaxis]).T * np.ones_like(CHI_VALS)
        PK_VALS = pk_interp_logk(CHI_VALS, np.log(K_VALS), grid=False)
        #compute integral via trapezium rule
        integrands = K1K2 * PK_VALS / CHI_VALS / CHI_VALS
        DCHI = CHI_VALS[:,1:] - CHI_VALS[:,:-1]
        integrands = DCHI * 0.5 * (integrands[:,1:] + integrands[:,:-1])
        c_ells = np.sum(integrands, axis=1)
        c_ell_errs = np.nan * np.ones_like(c_ells)

    elif method == "spline":
        for i,ell in enumerate(ells):
            nu = ell+0.5
            k_vals = nu / chi_vals
            pk_vals = pk_interp_logk(chi_vals, np.log(k_vals), grid=False)
            integrand = k1k2 * pk_vals / chi_vals / chi_vals
            int_spline = IUS(chi_vals, integrand)
            c_ells[i], c_ell_errs[i] = int_spline.integral(chimin, chimax), np.nan
    return c_ells, c_ell_errs
