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
         'external']


def ensure_starts_at_zero(z, nz, verbose):
    nbin = nz.shape[0]
    if z[0] > 0.00000001:
        z_new = np.zeros(len(z) + 1)
        z_new[1:] = z
        nz_new = np.zeros((nbin, len(z) + 1))
        nz_new[:, 1:] = nz
        nz_new[nz_new < 0] = 0
        if verbose:
            print("Putting n(0) = 0 at the start of the n(z)")
    else:
        z_new = z
        nz_new = nz
        nz_new[nz_new < 0] = 0

    return z_new, nz_new


def setup(options):
    print('\n')
    mode = options[option_section, 'mode']

    if mode not in modes:
        raise RuntimeError('Invalid mode set in config file. Please set one of\
        (inv-chi-)unified, (inv-chi-)separate, random or external')

    n_realisations = options[option_section, 'n_realisations']
    nz_filename = options[option_section, 'nz_file']
    data_set = options.get_string(option_section, "data_set")
    n_bins = options[option_section, 'n_bins']
    n_hist = options[option_section, 'n_hist']
    verbose = options.get_bool(option_section, "verbose", False)

    nz_file = pyfits.open(nz_filename)

    nz = np.zeros([n_realisations, n_bins, n_hist])
    gchi = np.zeros([n_realisations, n_bins, n_hist])
    chi = np.zeros([n_realisations, n_bins, n_hist])
    nz_mean = np.zeros([n_realisations, n_bins])
    inv_chi_mean = np.zeros([n_realisations, n_bins])

    for iext in np.arange(n_realisations):

        extname = 'nz_'+data_set+'_realisation_{0}'.format(iext)
        ext = nz_file[extname]
        zmid = ext.data['Z_MID']

        for ibin in np.arange(1, n_bins+1):
            nz[iext, ibin-1] = ext.data['BIN{0}'.format(ibin)]
            nz_mean[iext, ibin-1] = np.trapz(nz[iext, ibin-1]*zmid, zmid)
            if mode.startswith('inv-gchi'):
                chi, gchi[iext, ibin-1] = nz_to_gchi(zmid, nz[iext, ibin-1])
                inv_chi_mean[iext, ibin-1] = np.trapz(nz[iext, ibin-1]/chi,
                                                      chi)

    # Empty arrays for fiducial distribution data
    nz_fid = np.zeros([n_bins, n_hist])
    nz_mean_fid = np.zeros(n_bins)
    gchi_fid = np.zeros([n_bins, n_hist])
    chi_fid = np.zeros([n_bins, n_hist])
    inv_chi_mean_fid = np.zeros(n_bins)

    fiducial = True
    try:
        ext_fid = nz_file['NZ_SOURCE']
        zmid_fid = ext_fid.data['Z_MID']

    except:
        try:
            ext_fid = nz_file['nz_source']
            zmid_fid = ext_fid.data['Z_MID']
        except:
            if verbose:
                print('NZ_SOURCE fiducial extension not found. I will not \
                output expected position in the hyperrank')
                fiducial = False

    # If the fiducial extension exists, obtain all its properties.
    if fiducial:
        for ibin in np.arange(1, n_bins+1):
            nz_fid[ibin-1] = ext_fid.data['BIN{0}'.format(ibin)]
            nz_mean_fid[ibin-1] = np.trapz(nz_fid[ibin-1]*zmid_fid, zmid_fid)
            if mode.startswith('inv-gchi'):
                chi_fid, gchi_fid[ibin-1] = nz_to_gchi(zmid_fid, nz_fid[ibin-1])
                inv_chi_mean_fid[ibin-1] = np.trapz(nz_fid[ibin-1]/chi_fid, chi_fid)

    if mode == 'unified':
        fiducial_position = np.empty(1)
        nz_mean = nz_mean.mean(axis=1)
        nz_mean_fid = nz_mean_fid.mean()

        order = np.argsort(nz_mean)
        rank = np.argsort(order)

        ranked_nz = nz[order]
        ranked_nz_mean = nz_mean[order]

        if fiducial:
            fiducial_position[0] = np.searchsorted(ranked_nz_mean, nz_mean_fid)
            if verbose:
                print('Fiducial nz would have been ranked in the {0}th position out of {1}'.format(fiducial_position[0], n_realisations))

        if verbose:
            print('Ranked {0} unified distributions using mean redshift'.format(n_realisations))

    if mode == 'separate':
        fiducial_position = np.empty(n_bins)
        ranked_nz = np.empty([n_realisations, n_bins, n_hist])
        rank = np.empty([n_realisations, n_bins])
        for ibin in np.arange(1, n_bins+1):

            nz_mean_bin = nz_mean[:, ibin-1]

            order = np.argsort(nz_mean_bin)
            rank[:, ibin-1] = np.argsort(order)

            ranked_nz[:, ibin-1] = nz[order, ibin-1]
            ranked_nz_mean = nz_mean[order, ibin-1]

            if fiducial:
                fiducial_position[ibin-1] = np.searchsorted(ranked_nz_mean, nz_mean_fid[ibin-1])
                if verbose:
                    print('Fiducial bin {0} would have been ranked in the {1}th position out of {2}'.format(ibin, fiducial_position[ibin-1], n_realisations))

        if verbose:
            print('Ranked {0} separated tomographic bin distributions using mean redshift'.format(n_realisations))

    # External mode reads a list of ranks mapped to a index of redshift distributons stored in the input data set.
    if mode == 'external':
        external_filename = options[option_section, 'external_ranking_filename']
        if external_filename == "":
            raise RuntimeError('Set external mode but no input file defined. Aborting')

        external = np.genfromtxt(external_filename)
        order = np.argsort(external)
        rank = np.argsort(order)

        ranked_nz = nz[order]
        if verbose:
            print('Reading rankings from {0} '.format(external_filename))

        # Since we don't know which method was used to rank externally, we ignore this.
        fiducial = False

    if mode == 'random':
        # How do we keep the same ordering for all three processes? For now, just use the same random seed for all processes
        np.random.seed(n_realisations)

        order = np.arange(n_realisations)
        np.random.shuffle(order)
        rank = np.argsort(order)
        ranked_nz = nz[order]

        # No order, can't sort the fiducial distribution
        fiducial = False

    if mode == 'inv-gchi-unified':
        fiducial_position = np.empty(1)
        inv_chi_mean = inv_chi_mean.mean(axis=1)
        inv_chi_mean_fid = inv_chi_mean_fid.mean()

        order = np.argsort(inv_chi_mean)
        rank = np.argsort(order)

        ranked_nz = nz[order]
        ranked_nz_inv_chi_mean = inv_chi_mean[order]

        if fiducial:
            fiducial_position[0] = np.searchsorted(ranked_nz_inv_chi_mean, inv_chi_mean_fid)
            if verbose:
                print('Fiducial nz would have been ranked in the {0}th position'.format(fiducial_position[0]))

        if verbose:
            print('Ranked {0} unified distributions using mean inverse chi'.format(n_realisations))

    if mode == 'inv-gchi-separate':
        fiducial_position = np.empty(n_bins)
        ranked_nz = np.empty([n_realisations, n_bins, n_hist])
        rank = np.empty([n_realisations, n_bins])
        for ibin in np.arange(1, n_bins+1):

            inv_chi_mean_bin = inv_chi_mean[:, ibin-1]

            order = np.argsort(inv_chi_mean_bin)
            rank[:, ibin-1] = np.argsort(order)

            ranked_nz[:, ibin-1] = nz[order, ibin-1]
            ranked_nz_inv_chi_mean = inv_chi_mean[order, ibin-1]

            if fiducial:
                fiducial_position[ibin-1] = np.searchsorted(ranked_nz_inv_chi_mean, inv_chi_mean_fid[ibin-1])
                if verbose:
                    print('Fiducial nz would have been ranked in the {0}th position'.format(fiducial_position[ibin-1]))

        if verbose:
            print('Ranked {0} separated tomographic bin distributions using mean inverse chi'.format(n_realisations))

    config = {}
    config['sample'] = data_set.upper()
    config['ranked_nz'] = ranked_nz
    config['mode'] = mode
    config['n_realisations'] = n_realisations
    config['zmid'] = zmid
    config['n_bins'] = n_bins
    config['n_hist'] = n_hist
    config['verbose'] = verbose
    config['nsamples'] = 0
    config['fiducial'] = fiducial
    config['ranks'] = rank
    config['order'] = order
    if fiducial:
        config['fiducial_position'] = fiducial_position

    return config


def execute(block, config):

    ranked_nz = config['ranked_nz']
    mode = config['mode']
    pz = 'NZ_' + config['sample']
    n_realisations = config['n_realisations']
    verbose = config['verbose']
    nbin = config['n_bins']
    nhist = config['n_hist']
    block[pz, 'nbin'] = config['n_bins']
    ranks = config['ranks']

    if mode in ['unified', 'random', 'external', 'inv-gchi-unified']:

        sample = block['ranks', 'rank_hyperparm_1']
        sample = sample*n_realisations
        sample = int(sample)

        z, nz_sampled = ensure_starts_at_zero(config['zmid'], ranked_nz[sample], verbose)
        block[pz, 'z'] = z
        block[pz, 'nz'] = len(z)

        if verbose:
            print('Loaded ' + mode + ' ' + pz + '_{0} for rank {1}'.format(int(np.where(ranks == sample)[0]), sample))
            # print(orders)
            # print(ranks)

        for i in range(1, nbin+1):
            # write the nz to the data block. Is pz the name of the data_set
            block[pz, 'bin_{0}'.format(i)] = nz_sampled[i-1]

    if mode in ['separate', 'inv-gchi-separate']:
        # choose which ranked nz to use
        rank_indices = []
        nz_sampled = np.empty([nbin, nhist])
        for i in range(1, nbin+1):
            rank = block['ranks', 'rank_hyperparm_{0}'.format(i)]
            nz_sampled[i-1] = ranked_nz[int(rank*n_realisations), i-1]
            rank_indices.append(int(rank*n_realisations))

        z, nz_sampled = ensure_starts_at_zero(config['zmid'], nz_sampled, verbose)
        block[pz, 'z'] = z
        block[pz, 'nz'] = len(z)

        if verbose:
            print('Loaded separed bins from the following realisations:')
            for i in range(1, nbin+1):
                print('Bin {0}: '.format(i) + pz + '_{0}'.format(int(np.where(ranks[:, i-1] == rank_indices[i-1])[0])))

        for i in range(1, nbin+1):
            block[pz, 'bin_{0}'.format(i)] = nz_sampled[i-1]

    return 0


def cleanup(config):

    return 0
