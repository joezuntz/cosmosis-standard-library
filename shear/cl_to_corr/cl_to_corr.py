import scipy.interpolate
import pyfftlog
import numpy as np
from cosmosis.datablock import option_section

TRANSFORM_W = "w"
TRANSFORM_XI = "xi"
TRANSFORM_GAMMAT = "gamma"
TRANSFORM_XIP = "xip"
TRANSFORM_XIM = "xim"

TRANSFORMS = [TRANSFORM_W, TRANSFORM_GAMMAT, TRANSFORM_XI]

DEFAULT_SECTIONS = {
    TRANSFORM_XI:     ("shear_cl",        "shear_xi"),
    TRANSFORM_XIP:     ("shear_cl",        "shear_xi"),
    TRANSFORM_XIM:     ("shear_cl",        "shear_xi"),
    TRANSFORM_GAMMAT: ("galaxy_shear_cl", "galaxy_shear_xi"),
    TRANSFORM_W:      ("galaxy_cl",       "galaxy_xi"),
}
DEFAULT_N_TRANSFORM = 8192
DEFAULT_ELL_MIN = 0.0001
DEFAULT_ELL_MAX = 5.0e6
DEFAULT_THETA_MIN=0.1
DEFAULT_THETA_MAX=1000.0



# Bias q and order mu parameters for transform
_TRANSFORM_PARAMETERS = {
    TRANSFORM_W:       ( 0.0, 0.0),
    TRANSFORM_XIP:     ( 0.0, 0.0),
    TRANSFORM_XIM:     ( 0.0, 4.0),
    TRANSFORM_GAMMAT:  ( 0.0, 2.0),
}


class LogInterp(object):
    def __init__(self,angle,spec,kind):
        if np.all(spec>0):
            self.interp_func=scipy.interpolate.interp1d(np.log(angle),np.log(spec),kind,bounds_error=False)
            self.interp_type='loglog'
        elif np.all(spec<0):
            self.interp_func=scipy.interpolate.interp1d(np.log(angle),np.log(-spec),kind,bounds_error=False)
            self.interp_type='minus_loglog'
        else:
            self.interp_func=scipy.interpolate.interp1d(np.log(angle),spec,kind,bounds_error=False)
            self.interp_type="log_ang"

    def __call__(self,angle):
        if self.interp_type=='loglog':
            spec=np.exp(self.interp_func(np.log(angle)))
        elif self.interp_type=='minus_loglog':
            spec=-np.exp(self.interp_func(np.log(angle)))
        else:
            assert self.interp_type=="log_ang"
            spec=self.interp_func(np.log(angle))
        return spec



class Transformer(object):
    def __init__(self, transform_type, n, ell_min, ell_max, 
        theta_min, theta_max, lower=1.0, upper=-2.0):
        ell = np.logspace(np.log10(ell_min), np.log10(ell_max), n)
        dlogr = np.log(ell[1]) - np.log(ell[0])
        kropt = 1
        self.ell_min = ell_min
        self.ell_max = ell_max
        self.ell = ell
        self.q, self.mu = _TRANSFORM_PARAMETERS[transform_type]
        self.kr, self.xsave = pyfftlog.fhti(n, self.mu, dlogr, q=self.q, kropt=kropt)
        self.direction = -1
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.lower = lower
        self.upper = upper

        nc = 0.5*(n+1)
        log_ellmin = np.log(ell_min)
        log_ellmax = np.log(ell_max)
        log_ellmid = 0.5*(log_ellmin+log_ellmax)
        ell_mid = np.exp(log_ellmid)
        r_mid = self.kr/ell_mid  #radians
        x = np.arange(n)
        self.theta_rad = np.exp((x-nc)*dlogr) * r_mid #radians
        theta_arcmin = np.degrees(self.theta_rad) * 60.0 #arcmin
        self.range = (theta_arcmin>self.theta_min) & (theta_arcmin<self.theta_max)
        self.theta_arcmin = theta_arcmin[self.range]


    def __call__(self, ell_in, cl_in):
        """Convert the input ell and cl points to the points this transform requires, and then
        transform."""

        cl = self._interpolate_and_extrapolate_cl(ell_in, cl_in)

        if self.q==0:
            xi = pyfftlog.fht(self.ell*cl, self.xsave, tdir=self.direction)  / (2*np.pi) / self.theta_rad
        else:
            xi = pyfftlog.fhtq(self.ell*cl, self.xsave, tdir=self.direction)  / (2*np.pi) / self.theta_rad
        return self.theta_arcmin, xi[self.range]

    def _interpolate_and_extrapolate_cl(self, ell, cl):
        """Extrapolate and interpolate the input ell and cl to the default points for this transform"""
        ell_min = ell[0]
        ell_max = ell[-1]
        interpolator = LogInterp(ell, cl, 'linear')
        cl_out = interpolator(self.ell)
        bad_low = np.isnan(cl_out) & (self.ell<ell_min)
        bad_high = np.isnan(cl_out) & (self.ell>ell_max)
        
        cl_out[bad_low]  = cl[0]  * (self.ell[bad_low]  / ell_min)**self.lower
        cl_out[bad_high] = cl[-1] * (self.ell[bad_high] / ell_max)**self.upper
        
        return cl_out

class CosmosisTransformer(Transformer):
    def __init__(self, corr_type, options):
        # The type of transform to perform

        # Where to get/put the input/outputs
        default_input, default_output = DEFAULT_SECTIONS[corr_type]
        self.input_section = options.get_string(option_section, "input_section_name", default_input)
        self.output_section = options.get_string(option_section, "output_section_name", default_output)

        # Parameters of the transform
        n = options.get_int(option_section, "n_transform", DEFAULT_N_TRANSFORM)
        ell_min = options.get_double(option_section, "ell_min_extrapolate", DEFAULT_ELL_MIN)
        ell_max = options.get_double(option_section, "ell_max_extrapolate", DEFAULT_ELL_MAX)
        theta_min = options.get_double(option_section, "theta_min", DEFAULT_THETA_MIN)
        theta_max = options.get_double(option_section, "theta_max", DEFAULT_THETA_MAX)

        self.input_name = "bin_{}_{}"
        self.output_name = "bin_{}_{}"

        super(CosmosisTransformer, self).__init__(corr_type, n, ell_min, ell_max, theta_min, theta_max)

    def __call__(self, block):
        if block.has_value(self.input_section, "nbin_a"):
            nbin_a = block[self.input_section, "nbin_a"]
            nbin_b = block[self.input_section, "nbin_b"]
        else:
            nbin_a = nbin_b = block[self.input_section, "nbin"]
        ell = block[self.input_section, "ell"]
        for i in xrange(nbin_a):
            for j in xrange(nbin_b):
                b1 = i+1
                b2 = j+1
                input_name = self.input_name.format(b1,b2)
                output_name = self.output_name.format(b1,b2)

                if not block.has_value(self.input_section, input_name):
                    continue

                cl = block[self.input_section, input_name]
                theta, xi = super(CosmosisTransformer,self).__call__(ell, cl)

                #Cosmosis wants theta in radians
                theta = np.radians(theta/60.)

                block[self.output_section, "theta"] = theta
                block[self.output_section, output_name] = xi


class XiTransformer(object):
    def __init__(self, *args, **kwargs):
        self.xip = CosmosisTransformer(TRANSFORM_XIP, *args, **kwargs)
        self.xim = CosmosisTransformer(TRANSFORM_XIM, *args, **kwargs)
        self.xip.output_name = "xiplus_{}_{}"
        self.xim.output_name = "ximinus_{}_{}"

    def __call__(self, block):
        self.xip(block)
        self.xim(block)


def setup(options):
    corr_type = options.get_string(option_section, "corr_type")
    if corr_type not in TRANSFORMS:
        raise ValueError("Parameter transform in cl_to_corr must be one of {}".format(", ".join(TRANSFORMS)))

    if corr_type == TRANSFORM_XI:
        transformer = XiTransformer(options)
    else:
        # The transformer object, which stores all the 
        transformer = CosmosisTransformer(corr_type, options)

    return transformer


def execute(block, config):
    transformer = config
    transformer(block)
    return 0


