# Code by Donnacha Kirk
# Edited by Simon Samuroff 07/2015
# Edited by Jessie Muir 03/2023

from builtins import zip
from builtins import range
import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import scipy.interpolate
import warnings


def gaussian(z, mu, sigma):
    g = np.exp(-0.5 * (z - mu)**2 / sigma**2)
    dz = z[1] - z[0]
    gs = g.sum()
    if gs != 0:
        g /= g.sum() * dz
    return g


def delta(z, z0):
    x = np.zeros_like(z)
    # location nearest z0 in z
    idx = (np.abs(z - z0)).argmin()
    dz = z[1] - z[0]
    x[idx] = 1.0 / dz
    return x



def photometric_error(z, Nz, sigma_z, bias):
    nz = len(z)
    output = np.zeros((nz, nz))
    # If sigma==0 then all we have is an interpolation according to the bias value.
    # Using interpolation with kind=0 like this should ensure that using sigma=0 and sigma=very small give the same result.
    if sigma_z < 1e-3 and bias != 0:
        warnings.warn(
            "WARNING: Using a very small sigma_z can be problematic when also using a bias. Consider using the photoz_bias module separately.")
    if sigma_z == 0:
        zmax = z.max()
        for i in range(nz):
            p = delta(z, z[i] - bias)
            output[:, i] = p * Nz[i]

    else:
        for i in range(nz):
            # This doesn't work properly when you have a tiny sigma and a bias
            # because your n(z) can fall between the cracks.
            p = gaussian(z, z[i] - bias, sigma_z * (1 + z[i]))
            output[:, i] = p * Nz[i]
    return output


def find_bins(z, nz_true, nbin, zmin=None, zmax=None):
    """
    Find redshift bins with equal numbers of galaxies in them. 
    Inputs:
    z - redshift array
    nz_true - array containing full n(z) distribution
    nbin - number of bins to split that into
    zmin - lower bound where to start divisions; defaults to min 
           value in z array
    zmax - lower bound where to start divisions; defaults to max
           value in z array
    """
    if zmax is None:
        zmax = z.max()

    if zmin is not None:
        w =(z>zmin)*(z<zmax) # only sum counts in specified range
    else:
        zmin=0.
        w=1.


    nz_true = nz_true*w
    nz_true = nz_true / nz_true.sum() * nbin
    cum = np.cumsum(nz_true)
    bin_edges = [zmin]
    for i in range(1, nbin):
        edge = np.interp(1.0 * i, cum, z)
        bin_edges.append(edge)
    bin_edges.append(zmax)
    return np.array(bin_edges)


def compute_bin_nz(z_prob_matrix, z, edges, ngal):
    NI = []
    nbin = len(edges) - 1
    dz = z[1] - z[0]
   
    for low, high in zip(edges[:-1], edges[1:]):
        w = np.where((z > low) & (z < high))[0]
        # Sum over all possible ztrue
        # Equivalent to marginalising p(zphot|ztrue) wrt ztrue
        ni = z_prob_matrix[w, :].sum(axis=0)

        # Normalise the n(z) in each redshift bin to 1 over the redshift range
        # of the survey
        ni_integral = np.trapz(ni,z)
    
        ni *= 1.0 / ni_integral
        assert(len(ni) == len(z))
        NI.append(ni)

    return NI 

def get_ngals_presmooth(nz_true,z,edges,ngaltot):
    """
    Given bin divisions, and total number of galaxies in sample
    figure out how many galaxies are in each redshift bin. 

    N.b. Evaluation done before smoothing, with the idea being
    that you divide the sample in photo-z estimates, and then
    the smoothing translates photo-z selection to true n(z) 
    distributions. 
    """
    integrals = []
    for low, high in zip(edges[:-1], edges[1:]):
        w = (z>low)*(z<high)
        wnz = w*nz_true
        integrals.append(np.trapz(wnz,z))
    integrals = np.array(integrals)
    inttot = integrals.sum()
    ngal_fracs = integrals/inttot
    ngalvals = ngaltot*ngal_fracs
    return ngalvals


def smail_distribution(z, alpha, beta, z0):
    return (z**alpha) * np.exp(-(z / z0)**beta)


def compute_nz(alpha, beta, z0, z, nbin, sigma_z, ngal, bias, input_z_edges=None,force_equal=False):
    """
    alpha, beta, z0 - smail distribution variables
    z - array of redshift values
    nbin - number of bins to divide 
    sigma_z - photmetric redshift uncertainty at z=0
    ngal - total number density of galaxies in sample
    bias - shift in n(z) distribution + or - in redshift
    input_z_edges - optional; list of z values where to divide redshift bins
         Either expect nbin+1 values (for bin edges) or 2 values, which will 
         specify zmin and zmax for bin divisions (this zmin and max
         are distinct from the input options which control range of z
         values tabulated in array
    force_equal - if True, independentently of bin division, specifies 
           that all z bins will have equal number densities
    """
    # Set up Smail distribution of z vector as the distribution of true redshifts of the galaxies, n(ztrue)
    nz_true = smail_distribution(z, alpha, beta, z0)

    # Multiply that by a Gaussian to get the probability distribution of the measured photo-z for each true redshift p(zphot|ztrue)
    # This gives a 2D probability distribution
    z_prob_matrix = photometric_error(z, nz_true, sigma_z, bias)

    if input_z_edges is None:
        edges = find_bins(z, nz_true, nbin)
    elif (input_z_edges.size==2):
        # equal numbers in bins, with nominal min and max set by input_z_edges
        edges = find_bins(z, nz_true, nbin, input_z_edges[0], input_z_edges[1])
    else:
        edges = input_z_edges

    
    bin_nz = compute_bin_nz(z_prob_matrix, z, edges,ngal)
    if force_equal:
        ngals = ngal/nbin*np.ones(nbin)
    else:
        ngals = get_ngals_presmooth(nz_true,z,edges,ngal)
    return edges, bin_nz, ngals


def setup(options):
    # Determine how z values are tabulated
    dz = options.get_double(option_section, "dz", default=0.01)
    zmax = options.get_double(option_section, "zmax", default=4.0)
    zmin = options.get_double(option_section, "zmin", default=0.)

    # how many redshift bins are we dividing distribution into?
    nbin = options.get_int(option_section, "nbin")

    # specify that number densities for bins are equal
    force_equal = options.get_bool(option_section,"enforce_equal_numbers",default=False)

    # where to find input parameters and how to store in datablock
    in_section = options.get_string(option_section, "input_section", default=section_names.number_density_params)
    out_section = options.get_string(option_section, "output_section", default=section_names.wl_number_density)

    # can optionally pass in bin divisions manually
    if options.has_value(option_section, 'z_edges'):
        val = options[option_section,'z_edges']
        input_z_edges= np.atleast_1d(np.array(val, dtype=float))
        if (input_z_edges.size != 2) and (input_z_edges.size - 1 != nbin):
            raise ValueError("Incompatible nbin and zedges in smail module: nbin={} and z_edges.size={}".format(nbin,input_z_edges.size))
    else:
        input_z_edges = None
                            
    return (dz, zmin, zmax, nbin, in_section, out_section, input_z_edges,force_equal)

def execute(block, config):
    (dz, zmin, zmax, nbin, params, nz_section, input_z_edges, force_equal) = config
    alpha = block[params, "alpha"]
    beta = block[params, "beta"]
    z0 = block[params, "z0"]
    sigma_z = block[params, "sigz"]
    ngal = block[params, "ngal"]
    bias = block.get(params, "bias")

    # Compute the redshift vector
    z = np.arange(zmin, zmax + dz / 2, dz)

    # Run the main code for getting n(z) in bins
    edges, bins, ngals = compute_nz(alpha, beta, z0, z, nbin, sigma_z, ngal, bias, input_z_edges,force_equal)

    # Save the results
    block[nz_section, "nbin"] = nbin
    block[nz_section, "nz"] = len(z)
    block[nz_section, "z"] = z

    # Loop through the bins
    for i, bin in enumerate(bins):
        # The bin numbering starts at 1
        b = i + 1
        name = "BIN_%d" % b
        # Save the bin edges as parameters
        block[nz_section, "EDGE_%d" % b] = edges[i]
        # And save the bin n(z) as a column
        block[nz_section, name] = bin
        # also save Ngal for each bin
        block[nz_section,"NGAL_%d" %b] = ngals[i]
    # Also save the upper limit to the top bin
    block[nz_section, "EDGE_%d" % (nbin + 1)] = edges[-1]

    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
