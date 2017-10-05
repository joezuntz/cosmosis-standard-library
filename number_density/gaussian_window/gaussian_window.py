from builtins import range
from cosmosis.datablock import option_section
from cosmosis.datablock import names
import numpy as np
import scipy.integrate

DZ = 0.01
ZMAX = 3.0

# Option setup part.  Read options from in ifile.


def setup(options):
    z = options[option_section, "z"]
    sigma = options[option_section, "sigma"]
    section = options.get_string(
        option_section, "output_section", names.wl_number_density)
    # in case just one bin, convert to arrays:
    if np.isscalar(z):
        z = [z]
    if np.isscalar(sigma):
        sigma = [sigma]
    assert len(z) == len(
        sigma), "sigma and z had different sizes in the ini file section for gaussian_window"

    return (z, sigma, section)


def gaussian(x, mu, sigma):
    norm = np.sqrt(2 * np.pi) * sigma
    return np.exp(-0.5 * (x - mu)**2 / sigma**2) / norm


def execute(block, config):
    (mu, sigma, section) = config
    nbin = len(mu)
    # Set up z array.  Make sure ZMAX is inclusive
    z = np.arange(0, ZMAX + DZ, DZ)

    # Save Z array and count info
    block[section, "Z"] = z
    block[section, "NZ"] = len(z)
    block[section, "NBIN"] = nbin

    for i in range(1, nbin + 1):
        # generate simple gaussian window
        nz_bin = gaussian(z, mu[i - 1], sigma[i - 1])
        # the bin may not quite go to zero before we get to the
        # edges so normalize it
        norm = scipy.integrate.trapz(nz_bin, z)
        nz_bin /= norm
        # Save n(z)
        block[section, "BIN_%d" % i] = nz_bin
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
