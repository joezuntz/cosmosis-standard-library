from cosmosis.datablock import names, option_section
import numpy as np
import lens_kern
import scipy.integrate
from scipy.interpolate import RectBivariateSpline, interp1d

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters
distances = names.distances


def powerspec(k, z, rbs):
    return rbs.ev(k, z)


def clint(z, l, hspline, chispline, wa, wb, rbs):
    k = (l + 0.5) / chispline(z)
    tmp = wa(z) * wb(z) * powerspec(k, z, rbs) / chispline(z) ** 2 / hspline(z)
    return tmp


def cl_limber_x(z_chi, p_kz, l, k1, k2=None, xmin=0.0, xmax=13000.):
    """ calculate the cross-spectrum at multipole l between kernels k1 and k2 in the limber approximation.
    the distance integral is performed from conformal distance xmin to xmax (both in Mpc). """
    if k2 == None:
        k2 = k1

    def integrand(x):
        z = z_chi(x)
        return 1. / x ** 2 * k1.w_lxz(l, x, z) * k2.w_lxz(l, x, z) * rbs.ev((l + 0.5) / x, z)

    return scipy.integrate.quad(integrand, xmin, xmax, limit=100)[0]


def cl_limber_z(chi_z, hspline, rbs, l, k1, k2=None, zmin=0.0, zmax=1100.):
    """ calculate the cross-spectrum at multipole l between kernels k1 and k2 in the limber approximation.
    the distance integral is performed from redshifts zmin to zmax. """

    #  TODO check the H factor.

    if k2 == None:
        k2 = k1

    def integrand(z):

        x = chi_z(z)
        return 1. / x ** 2 / hspline(z) * k1.w_lxz(l, x, z) * k2.w_lxz(l, x, z) * rbs.ev((l + 0.5) / x, z)

    return scipy.integrate.quad(integrand, zmin, zmax, limit=100)[0]


def setup(options):

    llmin = options.get_double(option_section, "llmin", default=1.)
    llmax = options.get_double(option_section, "llmax", default=3.5)
    dlnl = options.get_double(option_section, "dlnl", default=.1)
    blockname = options.get_string(option_section, "matter_power", default="matter_power_nl")
    # blockname = options.get_string(option_section, "matter_power", default="matter_power_nl")
    # maybe the suffix for saving the spectra
    # zmax (or take it from CAMB)
    # maybe choose between kappa and others
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
    d_m *= h0
    h /= h0
    xlss *= h0
    # now in units of h^{-1} Mpc or the inverse
    chispline = interp1d(zdist, d_m)
    hspline = interp1d(zdist, h)

    lkern = lens_kern.kern(zdist, omega_m, h0, xlss)

    zmax = 6.
    zmin = 0.
    cl = np.zeros(np.size(lbins))
    for i, l in enumerate(lbins):
        cl[i] = scipy.integrate.quad(
            clint, zmin, zmax, args=(l, hspline, chispline, lkern.w_interp, lkern.w_interp, rbs))[0]
        print cl[i]
        cl[i] = cl_limber_z(chispline, hspline, rbs, l, lkern, k2=None, zmin=zmin, zmax=zmax)
        print cl[i]
        print ""
    section = "limber_spectra"
    block[section, "cl_kappa"] = cl
    block[section, "l_kappa"] = lbins

    return 0
