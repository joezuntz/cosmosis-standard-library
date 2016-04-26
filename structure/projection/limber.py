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
        ("status", ct.c_int),
]

LIMBER_STATUS_OK =  0
LIMBER_STATUS_ZERO =  1
LIMBER_STATUS_NEGATIVE =  2

c_gsl_spline = ct.c_void_p

dirname = os.path.split(__file__)[0]
lib = ct.cdll.LoadLibrary(os.path.join(dirname, "../../shear/spectra/interface.so"))
lib.get_named_w_spline.restype = c_gsl_spline
lib.get_named_w_spline.argtypes = [ct.c_size_t, ct.c_char_p, ct.c_int, c_dbl_array, ct.c_double, ct.c_void_p]

lib.get_named_nchi_spline.restype = c_gsl_spline
lib.get_named_nchi_spline.argtypes = [ct.c_size_t, ct.c_char_p, ct.c_int, c_dbl_array, ct.c_void_p, ct.c_void_p]

lib.cmb_wl_kappa_kernel.restype = c_gsl_spline
lib.cmb_wl_kappa_kernel.argtypes = [ct.c_double, ct.c_double, c_gsl_spline]

lib.get_kernel_peak.restype = ct.c_double
lib.get_kernel_peak.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int]

lib.limber_integral.restype = c_gsl_spline
lib.limber_integral.argtypes = [ct.POINTER(c_limber_config), ct.c_void_p, ct.c_void_p, ct.c_void_p]

lib.load_interpolator_chi.restype = ct.c_void_p
lib.load_interpolator_chi.argtypes = [ct.c_size_t, ct.c_void_p, ct.c_char_p, ct.c_char_p, ct.c_char_p, ct.c_char_p]

c_power_scaling_function = ct.CFUNCTYPE(ct.c_double, ct.c_double, 
    ct.c_double, ct.c_double, ct.c_voidp)

lib.load_interpolator_chi_function.restype = ct.c_void_p
lib.load_interpolator_chi_function.argtypes = [ct.c_size_t, ct.c_void_p, 
ct.c_char_p, ct.c_char_p, ct.c_char_p, ct.c_char_p, c_power_scaling_function, ct.c_void_p]


lib.interp_2d.restype = ct.c_double
lib.interp_2d.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]


lib.destroy_interp_2d.restype = None
lib.destroy_interp_2d.argtypes = [ct.c_void_p]



def get_named_nchi_spline(block, section, nbin, z, a_of_chi, chi_of_z):
    return GSLSpline(lib.get_named_nchi_spline(block._ptr, section, nbin, z, a_of_chi, chi_of_z))

def get_named_w_spline(block, section, bin, z, chi_max, a_of_chi):
    "Compute a shear kernel W(chi) spline"
    return GSLSpline(lib.get_named_w_spline(block._ptr, section, bin, z, chi_max, a_of_chi))

def get_cmb_kappa_spline(chi_max, chi_star, a_of_chi):
    "Compute the CMB WL kernel W_cmb(chi) spline"
    return GSLSpline(lib.cmb_wl_kappa_kernel(chi_max, chi_star, a_of_chi))

def evaluate_power(power, k, z):
    return lib.interp_2d(k,z,power)


def free_power(power):
    lib.destroy_interp_2d(power)



def load_power_chi(block, chi_of_z, section, k_name, z_name, p_name):
    "Load P(k,z) and convert z -> chi"
    r = lib.load_interpolator_chi(block._ptr, chi_of_z, section, k_name, z_name, p_name)
    if not r:
        raise ValueError("Could not load power spectrum from section {0} (k:{1} z:{2} p:{3})".format(section, k_name, z_name, p_name))
    return r

def load_power_chi_function(block, chi_of_z, section, k_name, z_name, p_name, function, args):
    "Load P(k,z) and convert z -> chi and scale P->f(k,z)*P"
    r = lib.load_interpolator_chi_function(block._ptr, chi_of_z, section, k_name, z_name, p_name, function, args)
    if not r:
        raise ValueError("Could not load scaled power spectrum from section {0} (k:{1} z:{2} p:{3})".format(section, k_name, z_name, p_name))
    return r

def get_kernel_peak(WX, WY, nchi=100):
    "Get chi of maximum of kernel"
    return lib.get_kernel_peak(WX, WY, nchi)

def limber(WX, WY, P, xlog, ylog, ell, prefactor):
    config = c_limber_config()
    config.xlog = xlog
    config.ylog = ylog
    config.n_ell = len(ell)
    config.ell = np.ctypeslib.as_ctypes(ell)
    config.prefactor = prefactor
    config.status = 0
    spline_ptr = lib.limber_integral(ct.byref(config), WX, WY, P)
    if config.status == LIMBER_STATUS_ZERO:
        ylog = False
    elif config.status == LIMBER_STATUS_NEGATIVE:
        raise ValueError("Negative value of the Limber integral.")
    spline = GSLSpline(spline_ptr, xlog=xlog, ylog=ylog)
    return spline
