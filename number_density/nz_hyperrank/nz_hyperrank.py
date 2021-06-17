"""
The Hyperrank module replaces 'load_nz_fits' and 'photoz_bias' to allow
marginalization of redshift systematics using a pool of distributions rather
than nuisance parameters.

On setup, realisations are read from an input file, ranked based on a
chatacteristic value of the whole distribution or their individual tomographic
bins and returned as an array.

Modes are:
'no-rank':
    No ranking. NZ array is filled on the same order the realisations are
    provided in the input fits file

'unified':
    Realisations are ranked according to the combined mean redshift across
    all tomographic bins and stored in the NZ array

'separate':
    Individual tomographic bins are ranked according to their mean
    redshift and stored on the NZ.

'external':
    Realisations are ranked according to the values on an external file.
    Values on the external file can be either the desired rank, or any
    other characteristic value describing the realisations.

'random':
    A random sequence of numbers is used to rank the distributions.
    This sequence then must remain constant during the pipeline.
    How do we keep the same ordering for all spawned processes?
    For now, just use the same random seed for all.

'inv-chi-unified':
    Similar to unified, but this time realisations are ranked according
    to the combined mean inverse comoving distance across all tomographic
    bins.

'inv-chi-separate':
    Similar to separete, but this time individual tomographic bins are
    ranked according to their mean inverse comoving distance.

On execute, the sampled n(z) is chosen based on the value of rank_hyperparm_i
defined by the sampler, which in turn is mapped to the rankings defined in
setup.
"""

try:
    from cosmosis.datablock import option_section
except:
    option_section = "options"
import numpy as np
from nz_gz import nz_to_gchi

try:
    import astropy.io.fits as pyfits
except ImportError:
    try:
        import pyfits
    except ImportError:
        raise RuntimeError("You need astropy installed to use the module \
        nz-hyperrank; try running: pip install astropy.")

modes = ['unified',
         'separate',
         'inv-chi-unified',
         'inv-chi-separate',
         'random',
         'external',
         'no-rank']

def load_histogram_form(ext, bin, upsampling):
    # Load the various z columns.
    # The cosmosis code is expecting something it can spline
    # so  we need to give it more points than this - we will
    # give it the intermediate z values (which just look like a step
    # function)
    zlow = ext.data['Z_LOW']
    zhigh = ext.data['Z_HIGH']

    if upsampling == 1:
        z = ext.data['Z_MID']
        nz = ext.data['BIN{0}'.format(bin)]

    else:
        z = np.linspace(0.0, zhigh[-1], len(zlow) * upsampling)
        sample_bin = np.digitize(z, zlow) - 1
        nz = ext.data['BIN{0}'.format(bin)][sample_bin]

    norm = np.trapz(nz, z)
    nz /= norm

    return z, nz

def ensure_starts_at_zero(z, nz):
    nbin = nz.shape[0]
    if z[0] > 0.00000001:
        z_new = np.zeros(len(z) + 1)
        z_new[1:] = z
        nz_new = np.zeros((nbin, len(z) + 1))
        nz_new[:, 1:] = nz
        nz_new[nz_new < 0] = 0
    else:
        z_new = z
        nz_new = nz
        nz_new[nz_new < 0] = 0

    return z_new, nz_new


def setup(options):
    # Returns a dictionary with an array of the realisations ranked according
    # the ranking mode.
    # The dictionary also includes additional information of the array and
    # selected modes.

    mode = options[option_section, 'mode']

    if mode not in modes:
        raise RuntimeError('Invalid mode set in config file. Please set one of\
        (inv-chi-)unified, (inv-chi-)separate, random or external')

    # Read configuration from inifile
    n_realisations = options[option_section, 'n_realisations']
    nz_filename = options[option_section, 'nz_file']
    data_set = options.get_string(option_section, "data_set")
    n_bins = options[option_section, 'n_bins']
    n_hist = options[option_section, 'n_hist']
    upsampling = options.get_int(option_section, 'upsampling', 1)

    nz_file = pyfits.open(nz_filename)

    # Initialize arrays for characteristic values from realisations
    nz = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    gchi = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    chi = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    nz_mean = np.zeros([n_realisations, n_bins])
    inv_chi_mean = np.zeros([n_realisations, n_bins])

    # Read all extensions from the input nz file and obtain their mean redshifts
    # and mean inverse comoving distances for all tomographic bins.
    for iext in np.arange(n_realisations):

        extname = 'nz_'+data_set+'_realisation_{0}'.format(iext)
        ext = nz_file[extname]

        for ibin in np.arange(1, n_bins+1):
            zmid, nz[iext, ibin-1] = load_histogram_form(ext, ibin, upsampling)
            nz_mean[iext, ibin-1] = np.trapz(nz[iext, ibin-1]*zmid, zmid)
            if mode.startswith('inv-chi'):
                chi, gchi[iext, ibin-1] = nz_to_gchi(zmid, nz[iext, ibin-1])
                inv_chi_mean[iext, ibin-1] = np.trapz(nz[iext, ibin-1]/chi, chi)


    # Depending on the mode selected, ranking begins here based on the nz_mean
    # or inv_chi_mean arrays.
    if mode == 'no-rank':

        ranked_nz = nz
        ranked_nz_mean = nz_mean


    if mode == 'unified':
        nz_mean = nz_mean.mean(axis=1)

        order = np.argsort(nz_mean)
        rank = np.argsort(order)

        ranked_nz = nz[order]
        ranked_nz_mean = nz_mean[order]


    if mode == 'separate':
        ranked_nz = np.empty([n_realisations, n_bins, n_hist*upsampling])
        rank = np.empty([n_realisations, n_bins])
        for ibin in np.arange(1, n_bins+1):

            nz_mean_bin = nz_mean[:, ibin-1]

            order = np.argsort(nz_mean_bin)
            rank[:, ibin-1] = np.argsort(order)

            ranked_nz[:, ibin-1] = nz[order, ibin-1]
            ranked_nz_mean = nz_mean[order, ibin-1]


    if mode == 'external':
        external_filename = options[option_section, 'external_ranking_filename']
        if external_filename == "":
            raise RuntimeError('Set external mode but no input file defined. Aborting')

        external = np.genfromtxt(external_filename)
        order = np.argsort(external)
        rank = np.argsort(order)

        ranked_nz = nz[order]


    if mode == 'random':
        np.random.seed(n_realisations)

        order = np.arange(n_realisations)
        np.random.shuffle(order)
        rank = np.argsort(order)
        ranked_nz = nz[order]


    if mode == 'inv-chi-unified':
        inv_chi_mean = inv_chi_mean.mean(axis=1)

        order = np.argsort(inv_chi_mean)
        rank = np.argsort(order)

        ranked_nz = nz[order]
        ranked_nz_inv_chi_mean = inv_chi_mean[order]


    if mode == 'inv-chi-separate':
        ranked_nz = np.empty([n_realisations, n_bins, n_hist*upsampling])
        rank = np.empty([n_realisations, n_bins])
        for ibin in np.arange(1, n_bins+1):

            inv_chi_mean_bin = inv_chi_mean[:, ibin-1]

            order = np.argsort(inv_chi_mean_bin)
            rank[:, ibin-1] = np.argsort(order)

            ranked_nz[:, ibin-1] = nz[order, ibin-1]
            ranked_nz_inv_chi_mean = inv_chi_mean[order, ibin-1]

    # create config dictionary with ranked rel and return it
    config = {}
    config['sample'] = data_set.upper()
    config['ranked_nz'] = ranked_nz
    config['mode'] = mode
    config['n_realisations'] = n_realisations
    config['zmid'] = zmid
    config['n_bins'] = n_bins
    config['n_hist'] = n_hist*upsampling

    return config


def execute(block, config):
    # Stores the sampled redshift distribution in the datablock by reading the
    # sampled rank_hyperparm_i values and mapping the ranked distributions to the
    # range of values rank_hyperparm_i can take, [0, 1)
    # There are two families of rank modes: Unified modes and separate modes.


    ranked_nz = config['ranked_nz']
    mode = config['mode']
    pz = 'NZ_' + config['sample']
    n_realisations = config['n_realisations']
    nbin = config['n_bins']
    nhist = config['n_hist']

    block[pz, 'nbin'] = config['n_bins']

    if mode in ['no-rank', 'unified', 'random', 'external', 'inv-chi-unified']:
        # A single rank_hyperparm_1 values is required, since all realisations
        # are considered a fixed collection of tomographic bins.
        sample = block['ranks', 'rank_hyperparm_1']
        sample = sample*n_realisations
        sample = int(sample)

        z, nz_sampled = ensure_starts_at_zero(config['zmid'], ranked_nz[sample])
        block[pz, 'z'] = z
        block[pz, 'nz'] = len(z)

        # write the nz to the data block. pz is the name of the data_set
        for i in range(1, nbin+1):
            block[pz, 'bin_{0}'.format(i)] = nz_sampled[i-1]

    if mode in ['separate', 'inv-chi-separate']:
        # A rank_hyperparm_i values for each tomographic bin is required
        # Realisations are not considered fixed groups of tomographic bins.
        # i.e., tomographic bin mixing is allowed.
        rank_indices = []
        nz_sampled = np.empty([nbin, nhist])
        for i in range(1, nbin+1):
            rank = block['ranks', 'rank_hyperparm_{0}'.format(i)]
            nz_sampled[i-1] = ranked_nz[int(rank*n_realisations), i-1]
            rank_indices.append(int(rank*n_realisations))

        z, nz_sampled = ensure_starts_at_zero(config['zmid'], nz_sampled)
        block[pz, 'z'] = z
        block[pz, 'nz'] = len(z)

        # write the nz to the data block. pz is the name of the data_set
        for i in range(1, nbin+1):
            block[pz, 'bin_{0}'.format(i)] = nz_sampled[i-1]

    return 0


def cleanup(config):

    return 0
