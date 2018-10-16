import os
import ctypes as ct
import numpy as np
from gsl_wrappers import GSLSpline, GSLSpline2d, BICUBIC, BILINEAR
import sys

c_dbl_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
c_int_array = np.ctypeslib.ndpointer(dtype=np.int, ndim=1, flags="C_CONTIGUOUS")

if sys.version_info[0] < 3:
    ascii_string = ct.c_char_p
else:
    class ascii_string(object):
        @classmethod
        def from_param(cls, value):
            if isinstance(value, bytes):
                return value
            else:
                return value.encode('ascii')


class c_limber_config(ct.Structure):
    _fields_ = [
        ("n_ell", ct.c_int),
        ("ell", ct.POINTER(ct.c_double)),
        ("prefactor", ct.c_double),
        ("status", ct.c_int),
        ("absolute_tolerance", ct.c_double),
        ("relative_tolerance", ct.c_double),
        ("K", ct.c_double),
]

class cl_config(ct.Structure):
    _fields_ = [
        ("n_ell", ct.c_int),
        ("ell", ct.POINTER(ct.c_int)),
        ("prefactor", ct.c_double),
        ("status", ct.c_int),
        ("delta_range_factor", ct.c_double),
        ("log_nu_range", ct.c_double),
        ("absolute_tolerance", ct.c_double),
        ("relative_tolerance", ct.c_double),
]

LIMBER_STATUS_OK =  0
LIMBER_STATUS_ZERO =  1
LIMBER_STATUS_NEGATIVE =  2
LIMBER_STATUS_ERROR =  3

c_gsl_spline = ct.c_void_p

dirname = os.path.split(__file__)[0]
lib = ct.cdll.LoadLibrary(os.path.join(dirname, "src/spec_tools.so"))
lib.get_named_w_spline.restype = c_gsl_spline
lib.get_named_w_spline.argtypes = [ct.c_size_t, ascii_string, ct.c_int, c_dbl_array, ct.c_double, ct.c_void_p]

lib.get_named_nchi_spline.restype = c_gsl_spline
lib.get_named_nchi_spline.argtypes = [ct.c_size_t, ascii_string, ct.c_int, c_dbl_array, ct.c_void_p, ct.c_void_p]

lib.get_reduced_kernel.restype = c_gsl_spline
lib.get_reduced_kernel.argtypes = [ct.c_void_p, ct.c_void_p]

lib.cmb_wl_kappa_kernel.restype = c_gsl_spline
lib.cmb_wl_kappa_kernel.argtypes = [ct.c_double, ct.c_double, c_gsl_spline, ct.c_double]

lib.get_kernel_peak.restype = ct.c_double
lib.get_kernel_peak.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_int]

lib.limber_integral.restype = ct.c_int
lib.limber_integral.argtypes = [ct.POINTER(c_limber_config), ct.c_void_p, ct.c_void_p, 
                                ct.c_void_p, ct.c_void_p, ct.c_int, c_dbl_array]

lib.load_interpolator_chi.restype = ct.c_void_p
lib.load_interpolator_chi.argtypes = [ct.c_size_t, ct.c_void_p, ascii_string, ascii_string, ascii_string, ascii_string]

c_power_scaling_function = ct.CFUNCTYPE(ct.c_double, ct.c_double, 
    ct.c_double, ct.c_double, ct.c_voidp)

lib.load_interpolator_chi_function.restype = ct.c_void_p
lib.load_interpolator_chi_function.argtypes = [ct.c_size_t, ct.c_void_p, 
ascii_string, ascii_string, ascii_string, ascii_string, c_power_scaling_function, ct.c_void_p]


lib.interp_2d.restype = ct.c_double
lib.interp_2d.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]

lib.destroy_interp_2d.restype = None
lib.destroy_interp_2d.argtypes = [ct.c_void_p]

lib.sigma_crit.restype = None
lib.sigma_crit.argtypes = [ct.c_void_p, ct.c_void_p, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]


def get_cmb_kappa_spline(chi_max, chi_star, a_of_chi, K=0.0):
    "Compute the CMB WL kernel W_cmb(chi) spline"
    return GSLSpline(lib.cmb_wl_kappa_kernel(chi_max, chi_star, a_of_chi, K))

def evaluate_power(power, k, z):
    return lib.interp_2d(k,z,power)

def free_power(power):
    try:
        lib.destroy_interp_2d(power)
    except ct.ArgumentError as e:
        power.__del__()

def load_power_z(block, section, k_name, z_name, p_name):
    z,k,p = block.get_grid(section, z_name, k_name, p_name)
    return GSLSpline2d(z, k, p.T, spline_type=BICUBIC)

def load_power_growth_chi(block, chi_of_z, section, k_name, z_name, p_name, k_growth=1.e-3):
    z,k,p = block.get_grid(section, z_name, k_name, p_name)
    chi = chi_of_z(z)
    growth_spline = growth_from_power(chi, k, p, k_growth)
    power_spline = GSLSpline2d(chi, np.log(k), p.T, spline_type=BICUBIC)
    return power_spline, growth_spline

def growth_from_power(chi, k, p, k_growth):
    "Get D(chi) from power spectrum"
    growth_ind=np.where(k>k_growth)[0][0]
    growth_array = np.sqrt(np.divide(p[:,growth_ind], p[0,growth_ind], 
                    out=np.zeros_like(p[:,growth_ind]), where=p[:,growth_ind]!=0.))
    return GSLSpline(chi, growth_array)

def get_reduced_kernel(orig_kernel, d_of_chi):
    return GSLSpline(lib.get_reduced_kernel(orig_kernel, d_of_chi))

def get_reduced_kernel(orig_kernel, d_of_chi):
    return GSLSpline(lib.get_reduced_kernel(orig_kernel, d_of_chi))

def get_named_nchi_spline(block, section, nbin, z, a_of_chi, chi_of_z):
    return GSLSpline(lib.get_named_nchi_spline(block._ptr, section, nbin, z, a_of_chi, chi_of_z))

def get_named_w_spline(block, section, bin, z, chi_max, a_of_chi):
    "Compute a shear kernel W(chi) spline"
    return GSLSpline(lib.get_named_w_spline(block._ptr, section, bin, z, chi_max, a_of_chi))

def load_power_chi(block, chi_of_z, section, k_name, z_name, p_name):
    "Load P(k,z) and convert z -> chi"
    r = lib.load_interpolator_chi(block._ptr, chi_of_z, section, k_name, z_name, p_name)
    if not r:
        raise ValueError("Could not load power spectrum from section {0} (k:{1} z:{2} p:{3})".format(section, k_name, z_name, p_name))
    return r

def get_kernel_peak(WX, WY, nchi=500):
    "Get chi of maximum of kernel"
    return lib.get_kernel_peak(WX, WY, nchi)

def get_sigma_crit(WX, WY):
    sigma_crit = ct.c_double(0.)
    chi_weighted = ct.c_double(0.)
    lib.sigma_crit(WX, WY, ct.byref(sigma_crit), ct.byref(chi_weighted))
    return sigma_crit.value, chi_weighted.value

def limber(WX, WY, P, D_chi, ell, prefactor, rel_tol=1.e-3, abs_tol=0., K=0.0 ):
    """D_chi is growth factor spline"""
    config = c_limber_config()
    config.n_ell = len(ell)
    config.ell = np.ctypeslib.as_ctypes(ell)
    config.prefactor = prefactor
    config.status = 0
    config.absolute_tolerance = abs_tol
    config.relative_tolerance = rel_tol
    config.K = K
    c_ell_out = np.empty(len(ell))

    #limber integral requries reduced kernels
    WX_red = get_reduced_kernel(WX, D_chi)
    WY_red = get_reduced_kernel(WY, D_chi)
    status = lib.limber_integral(ct.byref(config), WX_red, WY_red, P._ptr, 
                                 D_chi, 0, c_ell_out)
    if config.status  == LIMBER_STATUS_ERROR:
        raise ValueError("Error running Limber integral.  More info above.")
    return c_ell_out
