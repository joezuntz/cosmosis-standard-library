import os
import ctypes as ct
import numpy as np
from gsl_wrappers import GSLSpline


c_dbl_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

class c_limber_config(ct.Structure):
    _fields_ = [
        ("xlog", ct.c_bool),
        ("ylog", ct.c_bool),
        ("n_ell", ct.c_int),
        ("ell", ct.POINTER(ct.c_double)),
        ("prefactor", ct.c_double),
    ]

c_gsl_spline = ct.c_void_p

dirname = os.path.split(__file__)[0]
lib = ct.cdll.LoadLibrary(dirname + "/../../shear/spectra/interface.so")
lib.get_w_spline_named.restype = c_gsl_spline
lib.get_w_spline_named.argtypes = [ct.c_size_t, ct.c_char_p, ct.c_int, c_dbl_array, ct.c_double, ct.c_void_p]

lib.get_nchi_spline_named.restype = c_gsl_spline
lib.get_nchi_spline_named.argtypes = [ct.c_size_t, ct.c_char_p, ct.c_int, c_dbl_array, ct.c_void_p, ct.c_void_p]

lib.limber_integral.restype = c_gsl_spline
lib.limber_integral.argtypes = [ct.POINTER(c_limber_config), ct.c_void_p, ct.c_void_p, ct.c_void_p]

lib.load_interpolator_chi.restype = ct.c_void_p
lib.load_interpolator_chi.argtypes = [ct.c_size_t, ct.c_void_p, ct.c_char_p, ct.c_char_p, ct.c_char_p, ct.c_char_p]

def get_w_spline(block, nz_section, bin, z, chi_max, a_of_chi):
    "Compute a shear kernel W(chi) spline"
    return GSLSpline(lib.get_w_spline_named(block._ptr, nz_section, bin, z, chi_max, a_of_chi))

def load_power_chi(block, chi_of_z, section, k_name, z_name, p_name):
    "Load P(k,z) and convert z -> chi"
    return lib.load_interpolator_chi(block._ptr, chi_of_z, section, k_name, z_name, p_name)

def get_nchi_spline(block, nz_section, nbin, z, a_of_chi, chi_of_z):
    return GSLSpline(lib.get_nchi_spline_named(block._ptr, nz_section, nbin, z, a_of_chi, chi_of_z))

def limber(WX, WY, P, xlog, ylog, ell, prefactor):
    config = c_limber_config()
    config.xlog = xlog
    config.ylog = ylog
    config.n_ell = len(ell)
    config.ell = np.ctypeslib.as_ctypes(ell)
    config.prefactor = prefactor
    spline = GSLSpline(lib.limber_integral(ct.byref(config), WX, WY, P), xlog=xlog, ylog=ylog)
    return spline