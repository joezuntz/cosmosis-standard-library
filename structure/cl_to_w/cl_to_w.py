import scipy.special
import scipy.interpolate
from numpy import log, exp, cos, pi
from cosmosis.datablock import option_section
import numpy as np

def log_interp(x, y):
    s = scipy.interpolate.interp1d(log(x), log(y))
    x0 = x[0]
    y0 = y[0]
    x1 = x[-1]
    y1 = y[-1]
    def interpolator(xi):
        w1 = xi==0
        w2 = (xi>0) & (xi<=x0)
        w3 = xi>=x1
        w4 = ~ (w1 | w2 | w3)

        y = np.zeros_like(xi)
        y[w2] = y0 * (x0/xi[w2])
        y[w3] = y1 * (x1/xi[w3])**3
        y[w4] = exp(s(log(xi[w4])))
        return y
    return interpolator


def cl_to_w(ell, c_ell, theta):
    cl_interp = log_interp(ell, c_ell)
    ell_max = int(ell.max())
    ell_max_integral=200000
    ell_sample=np.arange(ell_max_integral)*1.0
    c_ell_sample = np.zeros(ell_max_integral)
    c_ell_sample = cl_interp(ell_sample)
#    for i,ell_i in enumerate(ell_sample):
#        c_ell_sample[i] = cl_interp(ell_i)

    f = (2*ell_sample+1)/(4*pi)
    w = np.zeros_like(theta)
    for i,t in enumerate(theta):
        p_ell, _ = scipy.special.lpn(ell_max_integral-1, cos(t))        
        w[i] = (f * p_ell * c_ell_sample).sum()
    return w


def setup(options):
    theta_min = options.get_double(option_section, "theta_min")
    theta_max = options.get_double(option_section, "theta_max")
    n_theta = options.get_int(option_section, "n_theta")
    theta = np.logspace(np.log10(theta_min), np.log10(theta_max), n_theta, endpoint=True)
    # a = theta[:-1]
    # b = theta[1:]
    # theta = 2./3. * (b**3-a**3)/(b**2-a**2)



    config = {}
    config['theta'] = theta

    return config

def execute(block, config):
    theta = config['theta']
    ell = block['galaxy_cl', 'ell']
    c_ell = block['galaxy_cl', 'bin_1_1']
    w = cl_to_w(ell, c_ell, np.radians(theta/60))
    block['galaxy_xi', 'theta'] = theta
    block['galaxy_xi', 'bin_1_1'] = w
    return 0
