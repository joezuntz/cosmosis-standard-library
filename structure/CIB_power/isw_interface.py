'''

Compute CIB cl, using Hall model, the single-spectra-energy-distribution (SSED) model of the  from Hall et. al. 2010 (arxiv:0912.4315)



'''


from cosmosis.datablock import names, option_section
import numpy as np
import isw_kernel
import lens_kern
import scipy.integrate
from scipy.interpolate import RectBivariateSpline, interp1d, InterpolatedUnivariateSpline

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.

cosmo = names.cosmological_parameters
distances = names.distances


def cl_limber_x(z_chi, p_kz, l, k1, k2=None, xmin=0.0, xmax=13000.):
    """ calculate the cross-spectrum at multipole l between kernels k1 and k2 in the limber approximation.
    the distance integral is performed from conformal distance xmin to xmax (both in Mpc). """
    if k2 == None:
        k2 = k1

    def integrand(x):
        z = z_chi(x)
        return 1. / x ** 2 * k1.w_lxz(l, x, z) * k2.w_lxz(l, x, z) * self.p_kz(l / x, z)

    return scipy.integrate.quad(integrand, xmin, xmax, limit=100)[0]


def cl_limber_z(chi_z, hspline, rbs, l, k1, k2=None, zmin=0.0, zmax=1100.):
    """ calculate the cross-spectrum at multipole l between kernels k1 and k2 in the limber approximation.
    the distance integral is performed from redshifts zmin to zmax. """

    #  TODO check the H factor.

    if k2 == None:
        k2 = k1

    def integrand(z):
        x = chi_z(z)
        # print 'hspline',hspline(z)
        return 1. / x ** 2 * hspline(z) * k1.w_lxz(l, x, z) * k2.w_lxz(l, x, z) * rbs.ev((l + 0.5) / x, z)

    return scipy.integrate.quad(integrand, zmin, zmax, limit=100)[0]


def powerspec(k, z, rbs):
    return rbs.ev(k, z)


def clint(z, l, hspline, chispline, wa, wb, rbs):

    k = (l + 0.5) / chispline(z)
    tmp = wa(z) * wb(z) * powerspec(k, z, rbs) / chispline(z) ** 2. / hspline(z)
    return tmp


def setup(options):

    # L BINS
    llmin = options.get_double(option_section, "llmin", default=1.)
    llmax = options.get_double(option_section, "llmax", default=3.)
    dlnl = options.get_double(option_section, "dlnl", default=.1)
    zmin = options.get_double(option_section, "zmin", default=1e-2)
    zmax = options.get_double(option_section, "zmax", default=10.)
    # Frequency of the CIB
    # What matter power spectrum to use, linear Halofit etc
    blockname = options.get_string(option_section, "matter_power", default="matter_power_lin")

    print 'llmin = ', llmin
    print 'llmax = ', llmax
    print 'dlnl = ', dlnl
    print 'matter_power = ', blockname
    print ' '

 #  the cib has a lot of parameter for the spectra ssed etc we can get them optionally here

    # maybe the suffix for saving the spectra
    # zmax (or take it from CAMB)
    # maybe choose between kappa and others
    lbins = np.arange(llmin, llmax, dlnl)
    lbins = 10. ** lbins
    return (lbins, blockname, zmin, zmax)


def execute(block, config):
    # Just a simple rename for clarity.
    lbins, blockname, zmin, zmax = config

    # lbins = np.array([
    #                820, 840, 860, 880])

    # LOAD POWER

    zpower = block[blockname, "z"]
    kpower = block[blockname, "k_h"]
    powerarray = block[blockname, "p_k"].reshape([np.size(zpower), np.size(kpower)]).T
    rbs = RectBivariateSpline(kpower, zpower, powerarray)

    # Cosmological parameters
    # =======================
    omega_m = block[cosmo, "omega_m"]
    h0 = block[cosmo, "h0"]
    # =======================

    # Growth
    # =======================
    # TODO Load and use just one of the 2 f or D
    d_z = block['growth_parameters', "d_z"]
    f_z = block['growth_parameters', "f_z"]
    z_growth = block['growth_parameters', "z"]
    dzispline = InterpolatedUnivariateSpline(z_growth, d_z * (1. + z_growth), k=5)
    fzispline = InterpolatedUnivariateSpline(z_growth, f_z, k=5)
    # =======================

    # Distances
    # =======================
    # reverse them so they are going in ascending order
    h = block[distances, "h"]
    tmp = h[::-1]
    h = tmp

    xlss = block[distances, "chistar"]

    zdist = block[distances, "z"]
    tmp = zdist[::-1]
    zdist = tmp

    d_m = block[distances, "d_m"]
    tmp = d_m[::-1]
    d_m = tmp
    # =======================

    # These have dimensions of Mpc; change to h^{-1} Mpc
    d_m *= h0
    h /= h0
    xlss *= h0
    # now in units of h^{-1} Mpc or the inverse
    chispline = interp1d(zdist, d_m)
    hspline = interp1d(zdist, h)
    # =======================
    # DEFINE KERNELs
    ikern = isw_kernel.kern(h0, omega_m, zdist, chispline, dzispline)
    # lkern = lens_kern.kern(zdist, omega_m, h0, xlss)
    cl = np.zeros(np.size(lbins))

    # COMPUTE CLs
    # =======================

    for i, l in enumerate(lbins):

        cl[i] = cl_limber_z(chispline, hspline, rbs, l, k1=ikern, k2=ikern, zmin=0., zmax=10.)
    # =======================
    # SAVE IN DATABLOCK

    section = "limber_spectra"
    block[section, "cl_cib"] = cl
    block[section, "l_cib"] = lbins

    return 0
