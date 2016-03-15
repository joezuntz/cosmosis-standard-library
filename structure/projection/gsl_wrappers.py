import os
import ctypes as ct
import numpy as np

def find_gsl():
    try:
        gsl_lib_dir = os.environ['GSL_LIB']
    except KeyError:
        raise ValueError("Please set the environment variable GSL_LIB to the directory containing libgsl.so or libgsl.dylib (mac)")

    libfile_so = os.path.join(gsl_lib_dir, "libgsl.so")
    libfile_dylib = os.path.join(gsl_lib_dir, "libgsl.dylib")
    if os.path.exists(libfile_so):
        libfile = libfile_so
        cblas_libfile = os.path.join(gsl_lib_dir, "libgslcblas.so")
    elif os.path.exists(libfile_dylib):
        libfile = libfile_dylib
        cblas_libfile = os.path.join(gsl_lib_dir, "libgslcblas.dylib")
    else:
        raise ValueError("Looked for GSL libs but could not find them: tried {0} and  {1}".format(libfile_so, libfile_dylib))
    return libfile, cblas_libfile

def load_gsl(libfile=None):
    if libfile is None:
        libfile, cblas_libfile = find_gsl()
    cblas = ct.CDLL(cblas_libfile, mode=ct.RTLD_GLOBAL)
    gsl = ct.CDLL(libfile, mode=ct.RTLD_GLOBAL)
    return gsl

#global gsl library
gsl = load_gsl()
c_dbl_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
gsl.gsl_spline_alloc.restype = ct.c_void_p
gsl.gsl_spline_alloc.argtypes = [ct.c_void_p, ct.c_size_t]
gsl.gsl_spline_init.restype = ct.c_int
gsl.gsl_spline_init.argtypes = [ct.c_void_p, c_dbl_array, c_dbl_array, ct.c_size_t]
gsl.gsl_spline_eval_e.restype = ct.c_int
gsl.gsl_spline_eval_e.argtypes = [ct.c_void_p, ct.c_double, ct.c_void_p, ct.POINTER(ct.c_double)]
gsl.gsl_spline_free.restype = None
gsl.gsl_spline_free.argtypes = [ct.c_void_p]




LINEAR = ct.c_void_p.in_dll(gsl, "gsl_interp_linear")
POLYNOMIAL = ct.c_void_p.in_dll(gsl, "gsl_interp_linear")
CSPLINE = ct.c_void_p.in_dll(gsl, "gsl_interp_cspline")
CSPLINE_PERIODIC = ct.c_void_p.in_dll(gsl, "gsl_interp_cspline_periodic")
AKIMA = ct.c_void_p.in_dll(gsl, "gsl_interp_akima")
AKIMA_PERIODIC = ct.c_void_p.in_dll(gsl, "gsl_interp_akima_periodic")



class GSLSpline(object):
    def __init__(self, x, y=None, spline_type=AKIMA, xlog=False, ylog=False):
        "Build a spline from either two numpy x and y arrays or a single void pointer"
        if y is None:
            self._ptr = x
        else:
            self._ptr = self._make_spline(x, y, spline_type)
        self.xlog = xlog
        self.ylog = ylog

    @staticmethod
    def _make_spline(x, y, spline_type):
        n = len(x)
        spline_ptr = gsl.gsl_spline_alloc(spline_type, n)
        gsl.gsl_spline_init(spline_ptr, x, y, n)
        return spline_ptr

    def __del__(self):
        if gsl is not None:
            gsl.gsl_spline_free(self._ptr)


    @property
    def _as_parameter_(self):
        return self._ptr

    def _eval(self, x):
        if self.xlog:
            x = np.log(x)
        y = ct.c_double(0.0)
        status = gsl.gsl_spline_eval_e(self._ptr, x, None, ct.byref(y))
        if status:
            raise Exception("GSL ERROR: {0}".format(status))
        y = y.value
        if self.ylog:
            y = np.exp(y)
        return y

    def __call__(self, x):
        """
        Evaluate the spline that this class points to.
        """
        if np.isscalar(x):
            return self._eval(x)
        else:
            x = np.array(x)
            return np.array([self._eval(xi) for xi in x])
