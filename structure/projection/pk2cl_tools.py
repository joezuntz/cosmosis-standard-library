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

def nearest_power_of_2(x):
    #Find nearest, greater power of 2 to x. 
    return 1<<(x-1).bit_length()

def get_dlogchi(dchi, chimax):
    #Since our chi values need to be log-spaced (for the fftlog), this dchi
    #will correspond to the last chi increment. Hence our dlogchi is given by
    #dlogchi = log(chi_N) - log(chi_{N-1}) = log(chi_N / chi_{N-1})
    #We can substitute chi_N = chimax and chi_{N-1} = chi_N - dchi to get
    #dlogchi = log(chimax / (chimax - dchi))
    dlogchi = np.log(chimax / (chimax-dchi))
    return dlogchi

def exact_integral(ells, kernel1_interp, kernel2_interp, pk0_interp_logk, growth_interp,
                    chimin, chimax, dlogchi, chi_pad_upper=1., chi_pad_lower=1., 
                    verbose=True):
    """full integral is \int_0^\inf k^{-1} dk P(k,0) I_1(k) I_2(k)
    where I_1(k) = \int_0^{\inf} k dr_1 F_1(r_1) r^{-0.5} D(r_1) J_{l+0.5}(kr_1),
    and F_1(r_1) is the radial kernel for tracer 1.
    We want to use a fftlog for the I_1(k) calculation, so write it in the form 
    I_1(k) = \int_0^{\inf} k dr_1 f_1(r_1) J_{mu}(kr_1).
    So f_1(r_1_) = F_1(r_1) D(r_1) r^{-0.5}
    We actually do the integral in log(k), so calculate \int_0^\inf dlogk P(k,0) I_1(k) I_2(k).
    """
    q=0
    assert chimin>0.
    log_chimin, log_chimax = np.log(chimin), np.log(chimax)
    if verbose:
        print("padding chi values by e^%.2f/%.2f at lower/upper ends"%(chi_pad_lower,chi_pad_upper))
    log_chimin_padded, log_chimax_padded = log_chimin-chi_pad_lower, log_chimax+chi_pad_upper
    nchi_orig = np.ceil((log_chimax-log_chimin)/dlogchi).astype(int)
    nchi = nearest_power_of_2(nchi_orig) #use nchi that is a power of 2 for fast fft.
    log_chi_vals = np.linspace(log_chimin_padded, log_chimax_padded, nchi)
    chi_vals = np.exp(log_chi_vals)
    if verbose:
        print("chimin padded, chimax padded, nchi padded:", chi_vals[0], chi_vals[-1], len(chi_vals))
    growth_vals = growth_interp(chi_vals)
    kernel1_vals = kernel1_interp(chi_vals)
    auto=False
    if kernel2_interp is kernel1_interp:
        kernel2_vals = kernel1_vals
    else:
        kernel2_vals = kernel2_interp(chi_vals)

    cell = np.zeros_like(ells)
    for i_ell, ell in enumerate(ells):

        f1_vals = kernel1_vals * growth_vals * np.power(chi_vals, -0.5)
        k_vals, I_1 = fft_log(chi_vals, f1_vals, q, ell+0.5)
        if auto:
            I_2 = I_1
        else:
            f2_vals = kernel2_vals * growth_vals * np.power(chi_vals, -0.5)
            _, I_2 = fft_log(chi_vals, f2_vals, q, ell+0.5)
        logk_vals = np.log(k_vals)
        pk_vals = pk0_interp_logk(logk_vals)
        #Now we can compute the full integral \int_0^{\inf} k dk P(k,0) I_1(k) I_2(k)
        #We are values logspaced in k, so calculate as \int_0^{inf} dlog(k) P(k,0) I_1(k) I_2(k)
        integrand_vals = pk_vals * I_1 * I_2
        #Spline and integrate the integrand.
        integrand_interp = IUS(logk_vals, integrand_vals)
        integral = integrand_interp.integral(logk_vals.min(), logk_vals.max())
        cell[i_ell] = integral
    return cell

def limber_integral(ells, kernel1, kernel2, pk_interp_logk, chimin, chimax, dchi,
    method="trapz", verbose=False):
    if verbose:
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
