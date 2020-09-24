try:
    from cosmosis.datablock import option_section
except:
    option_section = "options"
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
from nz_gz import nz_to_gchi

try:
    import astropy.io.fits as pyfits
except ImportError:
    try:
        import pyfits
    except ImportError:
        raise RuntimeError("You need astropy installed to use the module \
        nz-hyperrank; try running: pip install astropy.")

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

def gridmorph(x, shape, bounds_sigma=3., k_norm=2):
    ''' Function to map the point set x to a uniform grid
    `x` is an array of shape (M,N) giving locations of M
        points in N-dimensional space
    `shape` is the shape of the desired grid.  Should have
        N elements whose product equals M
    `bounds_sigma`: the input points will be initially
        linearly remapped to a unit interval in each
        dimension that spans mean+-`bounds_sigma`*std
        of the population.
    `k_norm`: the distance metric will be Euclidean
        distance raised to this power.

    Returns: an (M,N) array giving the coordinates of
        all the corresponding grid points.
    '''
    if len(shape)!=x.shape[1] or np.product(shape)!=x.shape[0]:
        raise ValueError('shape array does not match input points')
    ctr  = np.mean(x,axis=0)
    span = 2 * bounds_sigma * np.std(x,axis=0)

    xx = (x-ctr) / span + 0.5   # Rescale data
    xx = np.clip(xx,0.,1.)  # Clip into unit N-cube

    # Make the grid, rescaling to fill (0,1) in each dim
    uu = (np.indices(shape) + 0.5)
    # Flatten the grid, transpose to put point index first
    uu = uu.reshape(len(shape),-1).T
    uu /= np.array(shape)

    # Now calculate the distance metric between all pairs
    # of input points and grid points

    dx = xx[:,np.newaxis,:] - uu[np.newaxis,:,:]
    norm = np.sum(dx*dx,axis=2)  # Euclidean distance-squared
    if k_norm!=2.:
        norm = np.power(norm,k_norm/2.)

    # Now construct the optimized mapping
    ix, iu = linear_sum_assignment(norm)
    # The routine says the x indices will be returned
    # ordered.  So we can just use the 2nd set

    return uu[iu]

def factors(f, dim=2):
    # Returns the squarest? possible factors for a given number of realisations
    if dim == 1:
        return f
    p = np.zeros(dim, dtype=int)
    pow = float(dim)
    s = int(f**(1/pow))
    # s+2 will yield wrong results for small numbers.
    # should be s+1 but I don't know hot to deal with numeric errors in **1/3
    for i in range(1,s+2):
        if f % i == 0:
            p[0] = i
            p[1:] = factors(int(f/i), dim= dim-1)
    return p

def setup(options):

    def get_arr_ints(x):
        a = options[option_section, x]
        a = np.array(a, dtype=int)
        return a

    # Read configuration from inifile
    nz_filename = options[option_section, 'nz_file']
    data_set = options.get_string(option_section, "data_set")
    upsampling = options.get_int(option_section, 'upsampling', 1)
    mode = options.get_string(option_section, 'mode', 'mean')
    external_info = options.get_string(option_section, 'external_info', '')

    bin_ranks = get_arr_ints('bin_ranks')
    dimensions = options.get_int(option_section, 'dimensions', 2)
    assert len(bin_ranks) == dimensions, 'You asked for {} dimensions but provided {} bins for the ranking.'.format(dimensions, len(bin_ranks))

    resume = options.get_bool(option_section, 'resume', False)
    resume_map = options.get_string(option_section, 'resume_map', "")

    # Determine number of realisations and array shapes from the file
    nz_file = pyfits.open(nz_filename)
    n_realisations = 0
    n_bins = 0
    istart = 0 # This is so we can index extensions by number instead of name. Much faster

    for iext in np.arange(1, len(nz_file)):
        if nz_file[iext].header['EXTNAME'].startswith('nz_{0}_realisation'.format(data_set)):
            n_realisations += 1
        if nz_file[iext].header['EXTNAME'].startswith('nz_{0}_realisation_0'.format(data_set)):
            istart = iext
            n_hist = len(nz_file[iext].data['Z_MID'])
            for col in nz_file[iext].data.columns:
                if col.name.startswith('BIN'):
                    n_bins += 1

    print('Multirank detected {0} realisations, {1} tomographic bins, {2} histogram bins.'.format(n_realisations, n_bins, n_hist))

    # Initialize arrays for characteristic values from realisations
    nz = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    gchi = np.zeros([n_realisations, n_bins, n_hist*upsampling])
    nz_mean = np.zeros([n_realisations, n_bins])
    inv_chi_mean = np.zeros([n_realisations, n_bins])

    for iext in np.arange(n_realisations):

        ext = nz_file[iext + istart]

        for ibin in np.arange(n_bins):
            zmid, nz[iext, ibin] = load_histogram_form(ext, ibin+1, upsampling)
            nz_mean[iext, ibin] = np.trapz(nz[iext, ibin]*zmid, zmid)
            if mode == 'invchi':
                chi, gchi[iext, ibin] = nz_to_gchi(zmid, nz[iext, ibin])
                inv_chi_mean[iext, ibin] = np.trapz(nz[iext, ibin]/chi, chi)

    if mode == 'mean':
        xx = nz_mean[:,bin_ranks-1]
    if mode == 'invchi':
        xx = inv_chi_mean[:,bin_ranks-1]
    if mode == 'external':
        xx = np.load(external_info)

    map_shape = factors(n_realisations, dim=dimensions)

    # Read previously computed uniform map or compute a new one.
    if resume:
        try:
            uu = np.load(resume_map)
            uu_loaded_dim = np.zeros(uu.shape[-1])
            for idim in range(uu.shape[-1]):
                uu_loaded_dim[idim] = len(np.unique(uu[:,idim]))
            print('Found uniform map of dimensions ' + str(uu_loaded_dim))
        except:
            print('Tried to resume using previous uniform map but could not find it.')
            print('Generating a new one of dimensions ' + str(map_shape) + ' with tomographic bins ' + str(bin_ranks))
            uu = gridmorph(xx, map_shape)
            np.save(resume_map, uu)
    else:
        print('Multirank is generating a {}-dimensional uniform map to sample realisations.'.format(dimensions))
        print('Dimensions are ' + str(map_shape) + ' and uses tomographic bins ' + str(bin_ranks))
        uu = gridmorph(xx, map_shape)

    assert dimensions == uu.shape[-1], "Loaded uniform map was generated with a different dimensionality."
    assert map_shape[0] == len(np.unique(uu[:,0])), "Loaded uniform map has a different shape than the one computed here."
    assert map_shape[1] == len(np.unique(uu[:,1])), "Loaded uniform map has a different shape than the one computed here."

    # create config dictionary with mapped realisations and return it
    config = {}
    config['sample'] = data_set.upper()
    config['nz'] = nz
    config['nz_mean'] = nz_mean
    config['zmid'] = zmid
    config['uniform_grid'] = uu
    config['dimensions'] = dimensions
    config['map_shape'] = map_shape

    return config

def execute(block, config):
    # Stores the sampled redshift distribution in the datablock by reading the
    # sampled rank_hyperparm_i values and mapping the ranked distributions to the
    # range of values rank_hyperparm_i can take, [0, 1)

    nz = config['nz']

    zmid = config['zmid']
    uu = config['uniform_grid']
    pz = 'NZ_' + config['sample']
    nz_mean = config['nz_mean']
    nbins = nz.shape[1]
    dimensions = config['dimensions']
    map_shape = config['map_shape']

    # Extract coordinates from sampled hyperparameters and find index from map

    ranks = [block['ranks', 'rank_hyperparm_{}'.format(idim+1)] for idim in range(dimensions)]
    coordinates = [(int(ranks[idim]*map_shape[idim]) + 0.5)/map_shape[idim] for idim in range(dimensions)]
    uu_index = [np.argwhere(uu[:,idim] == coordinates[idim]).flatten() for idim in range(dimensions)]

    index = uu_index[-1]
    for idim in range(dimensions-1):
        index = np.intersect1d(index, uu_index[idim], assume_unique=True)
    index = int(index)

    # Extract sampled nz
    z, nz_sampled = ensure_starts_at_zero(zmid, nz[index])

    # Write info to datablock
    block[pz, 'nbin'] = nbins
    block[pz, 'z'] = z
    block[pz, 'nz'] = len(z)

    for ibin in np.arange(nbins):
        block[pz, 'bin_{0}'.format(ibin+1)] = nz_sampled[ibin]
        block['ranks', 'mean_z_{0}'.format(ibin)] = nz_mean[index, ibin]

    block['ranks', 'realisation_id'.format(ibin)] = index

    return 0

def cleanup(config):

    return 0
