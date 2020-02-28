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
    mode = options[option_section, 'mode']

    if mode not in modes:
        raise RuntimeError('Invalid mode set in config file. Please set one of\
        (inv-chi-)unified, (inv-chi-)separate, random or external')

    n_realisations = options[option_section, 'n_realisations']
    nz_filename = options[option_section, 'nz_file']
    data_set = options.get_string(option_section, "data_set")
    n_bins = options[option_section, 'n_bins']
    n_hist = options[option_section, 'n_hist']
    upsampling = options.get_int(option_section, 'upsampling', 1)

    nz_file = pyfits.open(nz_filename)

    nz = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    gchi = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    chi = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    nz_mean = np.zeros([n_realisations, n_bins])
    inv_chi_mean = np.zeros([n_realisations, n_bins])

    for iext in np.arange(n_realisations):

        extname = 'nz_'+data_set+'_realisation_{0}'.format(iext)
        ext = nz_file[extname]

        for ibin in np.arange(1, n_bins+1):
            zmid, nz[iext, ibin-1] = load_histogram_form(ext, ibin, upsampling)
            nz_mean[iext, ibin-1] = np.trapz(nz[iext, ibin-1]*zmid, zmid)
            if mode.startswith('inv-chi'):
                chi, gchi[iext, ibin-1] = nz_to_gchi(zmid, nz[iext, ibin-1])
                inv_chi_mean[iext, ibin-1] = np.trapz(nz[iext, ibin-1]/chi, chi)


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
        # How do we keep the same ordering for all three processes? For now, just use the same random seed for all processes
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

    ranked_nz = config['ranked_nz']
    mode = config['mode']
    pz = 'NZ_' + config['sample']
    n_realisations = config['n_realisations']
    nbin = config['n_bins']
    nhist = config['n_hist']

    block[pz, 'nbin'] = config['n_bins']

    if mode in ['no-rank', 'unified', 'random', 'external', 'inv-chi-unified']:

        sample = block['ranks', 'rank_hyperparm_1']
        sample = sample*n_realisations
        sample = int(sample)

        z, nz_sampled = ensure_starts_at_zero(config['zmid'], ranked_nz[sample])
        block[pz, 'z'] = z
        block[pz, 'nz'] = len(z)

        for i in range(1, nbin+1):
            # write the nz to the data block. Is pz the name of the data_set
            block[pz, 'bin_{0}'.format(i)] = nz_sampled[i-1]

    if mode in ['separate', 'inv-chi-separate']:
        # choose which ranked nz to use
        rank_indices = []
        nz_sampled = np.empty([nbin, nhist])
        for i in range(1, nbin+1):
            rank = block['ranks', 'rank_hyperparm_{0}'.format(i)]
            nz_sampled[i-1] = ranked_nz[int(rank*n_realisations), i-1]
            rank_indices.append(int(rank*n_realisations))

        z, nz_sampled = ensure_starts_at_zero(config['zmid'], nz_sampled)
        block[pz, 'z'] = z
        block[pz, 'nz'] = len(z)

        for i in range(1, nbin+1):
            block[pz, 'bin_{0}'.format(i)] = nz_sampled[i-1]

    return 0


def cleanup(config):

    return 0
