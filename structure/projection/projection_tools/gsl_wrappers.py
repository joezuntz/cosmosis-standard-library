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
        try:
            libfile, cblas_libfile = find_gsl()
        except ValueError:
            libfile = "libgsl.so"
            cblas_libfile = "libgslcblas.so"
    # This isn't always necessary, but it is on NERSC
    # To help debug I'll leave this lib1 variable here.
    try:
        lib1 = ct.CDLL(cblas_libfile, mode=ct.RTLD_GLOBAL)
    except OSError:
        lib1 = None
    try:
        # Try loading the CBLAS file first. Sometimes this is needed.
        gsl = ct.CDLL(libfile, mode=ct.RTLD_GLOBAL)
    except OSError as err:
        msg = f"{err} ** If you are using conda you might need to source cosmosis-configure, or otherwise you may need to set up GSL differently. **"
        raise OSError(msg) from err
    return gsl

#global gsl library
gsl = load_gsl()
c_dbl_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

#1d splines
gsl.gsl_spline_alloc.restype = ct.c_void_p
gsl.gsl_spline_alloc.argtypes = [ct.c_void_p, ct.c_size_t]
gsl.gsl_spline_init.restype = ct.c_int
gsl.gsl_spline_init.argtypes = [ct.c_void_p, c_dbl_array, c_dbl_array, ct.c_size_t]
gsl.gsl_spline_eval_e.restype = ct.c_int
gsl.gsl_spline_eval_e.argtypes = [ct.c_void_p, ct.c_double, ct.c_void_p, ct.POINTER(ct.c_double)]
gsl.gsl_spline_free.restype = None
gsl.gsl_spline_free.argtypes = [ct.c_void_p]

#2d splines
gsl.gsl_spline2d_alloc.restype = ct.c_void_p
gsl.gsl_spline2d_alloc.argtypes = [ct.c_void_p, ct.c_size_t, ct.c_size_t]
gsl.gsl_spline2d_init.restype = ct.c_int
gsl.gsl_spline2d_init.argtypes = [ct.c_void_p, c_dbl_array, c_dbl_array, c_dbl_array, ct.c_size_t, ct.c_size_t]
gsl.gsl_spline2d_eval_e.restype = ct.c_int
gsl.gsl_spline2d_eval_e.argtypes = [ct.c_void_p, ct.c_double, ct.c_double, ct.c_void_p, ct.POINTER(ct.c_double)]
gsl.gsl_spline2d_free.restype = None
gsl.gsl_spline2d_free.argtypes = [ct.c_void_p]

LINEAR = ct.c_void_p.in_dll(gsl, "gsl_interp_linear")
POLYNOMIAL = ct.c_void_p.in_dll(gsl, "gsl_interp_linear")
CSPLINE = ct.c_void_p.in_dll(gsl, "gsl_interp_cspline")
CSPLINE_PERIODIC = ct.c_void_p.in_dll(gsl, "gsl_interp_cspline_periodic")
AKIMA = ct.c_void_p.in_dll(gsl, "gsl_interp_akima")
AKIMA_PERIODIC = ct.c_void_p.in_dll(gsl, "gsl_interp_akima_periodic")
BILINEAR = ct.c_void_p.in_dll(gsl, "gsl_interp2d_bilinear")
BICUBIC = ct.c_void_p.in_dll(gsl, "gsl_interp2d_bicubic")

class NullSplineError(ValueError):
    pass            

class GSLSpline(object):
    def __init__(self, x, y=None, spline_type=AKIMA, xlog=False, ylog=False):
        "Build a spline from either two numpy x and y arrays or a single void pointer"
        if y is None:
            self._ptr = x
            if x is None:
                raise NullSplineError("Tried to wrap a null pointer in GSLSpline")
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
        if gsl is not None and self._ptr is not None:
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

class GSLSpline2d(object):
    def __init__(self, x, y=None, Z=None, spline_type=BILINEAR):
        if (y is None) and (Z is None):
            self._ptr=x
            if x is None:
                raise NullSplineError("Tried to wrap a null pointer in GSLSpline")
        else:
            self._ptr = self._make_spline(x, y, Z, spline_type=spline_type)

    @staticmethod
    def _make_spline(x, y, Z, spline_type='bilinear'):
        nx = len(x)
        ny = len(y)
        #print x,y,Z
        #print 'allocating %d,%d spline'%(nx,ny)
        spline_ptr = gsl.gsl_spline2d_alloc(spline_type, nx, ny)
        #print 'initialising 2d spline'
        gsl.gsl_spline2d_init(spline_ptr, x, y, Z.flatten(), nx, ny)
        #print 'spline initialised'
        return spline_ptr

    def __del__(self):
        if gsl is not None and self._ptr is not None:
            gsl.gsl_spline2d_free(self._ptr)

    def _eval(self, x, y):
        z = ct.c_double(0.0)
        status = gsl.gsl_spline2d_eval_e(self._ptr, x, y, None, None, ct.byref(z))
        if status:
            raise Exception("GSL ERROR: {0}".format(status))
        z = z.value
        return z

    def __call__(self, x, y):
        """
        Evaluate the spline that this class points to - for testing, you wouldn't actually
        want to do this at the python level
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        if len(x)==1:
            return np.array([self._eval(x[0],yi) for yi in y])
        elif len(y)==1:
            return np.array([self._eval(xi, y[0]) for xi in x])
        else:
            return np.array([self._eval(xi,yi) for (xi,yi) in zip(x,y)])


def test1d():
    x = np.linspace(0,1,10)
    y = x**0.5
    spline = GSLSpline(x, y)
    output = ct.c_double(0.0)
    x_test = np.random.random()
    gsl.gsl_spline_eval_e(spline._ptr, x_test, None, ct.byref(output))
    print(x_test,output.value/x_test**0.5)

def test2d():
    x = np.linspace(0,1,10)
    y = np.linspace(0,1,10)
    X,Y = np.meshgrid(x,y)
    Z = X * Y**0.5
    spline = GSLSpline2d(x, y, Z)
    output = ct.c_double(0.0)
    x_test = np.random.random()
    y_test = np.random.random()
    print('calling 2d spline')
    gsl.gsl_spline2d_eval_e(spline._ptr, x_test, y_test, None, None, ct.byref(output))
    print(output.value, x_test*y_test**0.5)

if __name__=="__main__":
    test2d()
