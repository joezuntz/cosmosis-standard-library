from cosmosis.datablock import names, option_section
import numpy as np
import lens_kern
import scipy.integrate
from scipy.interpolate import RectBivariateSpline, interp1d, UnivariateSpline
# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters
distances = names.distances

def powerspec(k, z, rbs):
    return rbs.ev(k, z)


def clint(z, l, hspline, chispline, wa, wb, rbs):
    k = l / chispline(z)
    tmp = wa(z) * wb(z) * powerspec(k, z, rbs)  / chispline(z) ** 2 / hspline(z)
    return tmp


def setup(options):

    llmin = options.get_double(option_section, "llmin", default=1.)
    llmax = options.get_double(option_section, "llmax", default=3.5)
    dlnl = options.get_double(option_section, "dlnl", default=.1)
    print 'setup point'
    blockname = options.get_string(option_section, "matter_power", default="matter_power_lin")
    lbins = np.arange(llmin, llmax, dlnl)
    lbins = 10. ** lbins
    return (lbins, blockname)


def execute(block, config):
    # Just a simple rename for clarity.
    lbins, blockname = config

    # LOAD POWER

    zpower = block[blockname, "z"]
    kpower = block[blockname, "k_h"]
    powerarray = block[blockname, "p_k"].reshape([np.size(zpower), np.size(kpower)]).T
    rbs = RectBivariateSpline(kpower, zpower, powerarray)

    # Cosmological parameters

    omega_m = block[cosmo, "omega_m"]
    h0 = block[cosmo, "h0"]
    # Distances
    h = block[distances, "h"]
    tmp = h[::-1]
    h = tmp

    zdist = block[distances, "z"]
    tmp = zdist[::-1]  # reverse them so they are going in ascending order
    zdist = tmp

    d_m = block[distances, "d_m"]
    tmp = d_m[::-1]
    d_m = tmp

    xlss = block[distances, "chistar"]

    # These have dimensions of Mpc; change to h^{-1} Mpc
    d_m = d_m * h0
    h = h / h0
    xlss*=h0
    # now in units of h^{-1} Mpc or the inverse
    chispline = interp1d(zdist, d_m)
    hspline = interp1d(zdist, h)

    lkern = lens_kern.kern(zdist, omega_m, h0, xlss)

    zmax = 30.
    zmin = 0.
    cl = np.zeros(np.size(lbins))
    for i, l in enumerate(lbins):
        cl[i] = scipy.integrate.quad(clint, zmin, zmax, args=(l, hspline, chispline, lkern.w_interp, lkern.w_interp, rbs))[0]
    section = "IA_cls"
    block[section, "cl"] = cl
    block[section, "l"] = lbins


# Now go with the integral and save


    return 0
