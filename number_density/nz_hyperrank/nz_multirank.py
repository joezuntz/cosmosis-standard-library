try:
    from cosmosis.datablock import option_section
except:
    option_section = "options"
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt

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

def setup(options):

    # Read configuration from inifile
    nz_filename = options[option_section, 'nz_file']
    data_set = options.get_string(option_section, "data_set")
    upsampling = options.get_int(option_section, 'upsampling', 1)
    debug = options.get_bool(option_section, 'debug', False)

    dim1 = options.get_int(option_section, 'dimension_1')
    dim2 = options.get_int(option_section, 'dimension_2')
    bin_rank1 = options.get_int(option_section, 'bin_rank1', 1)
    bin_rank2 = options.get_int(option_section, 'bin_rank2', 4)

    nz_file = pyfits.open(nz_filename)

    n_realisations = 0
    n_bins = 0
    istart = 0

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
    nz_mean = np.zeros([n_realisations, n_bins])

    for iext in np.arange(n_realisations):

        ext = nz_file[iext + istart]

        for ibin in np.arange(n_bins):
            zmid, nz[iext, ibin] = load_histogram_form(ext, ibin+1, upsampling)
            nz_mean[iext, ibin] = np.trapz(nz[iext, ibin]*zmid, zmid)

    xx = nz_mean[:,[bin_rank1-1, bin_rank2-1]]
    if debug:
        print('Multirank is generating a uniform map to sample your realisations')
    uu = gridmorph(xx, [dim1, dim2])

    if debug:
        plt.plot(xx[:,0],xx[:,1],'ro')
        plt.gca().set_aspect('equal')
        plt.title('original point positions')
        plt.savefig('original_positions.png', dpi=200)

        f, ax = plt.subplots(1,2, figsize=(10,4))

        ax[0].scatter(uu[:,0],uu[:,1],c=xx[:,0],cmap='Spectral', s=5)
        ax[0].set_aspect('equal')
        ax[0].set_title('Original mean z {}'.format(bin_rank1))

        ax[1].scatter(uu[:,0],uu[:,1],c=xx[:,1],cmap='Spectral', s=5)
        ax[1].set_aspect('equal')
        ax[1].set_title('Original z {}'.format(bin_rank2))

        plt.savefig('uniform_grid_means.png', dpi=200)

    # create config dictionary with ranked realisations and return it
    config = {}
    config['sample'] = data_set.upper()
    config['nz'] = nz
    config['nz_mean'] = nz_mean
    config['zmid'] = zmid
    config['uniform_grid'] = uu
    config['dimensions'] = (dim1, dim2)
    config['debug'] =  debug

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
    dim1, dim2 = config['dimensions']
    debug = config['debug']

    # Extract coordinates from sampled hyperparameters

    rank_1 = block['ranks', 'rank_hyperparm_1']
    rank_2 = block['ranks', 'rank_hyperparm_2']
    coordinate_1 = (int(rank_1*dim1) + 0.5)/dim1
    coordinate_2 = (int(rank_2*dim2) + 0.5)/dim2
    uu_1 = np.argwhere(uu[:,0] == coordinate_1).flatten()
    uu_2 = np.argwhere(uu[:,1] == coordinate_2).flatten()

    index = int(np.intersect1d(uu_1, uu_2, assume_unique=True))
    z, nz_sampled = ensure_starts_at_zero(zmid, nz[index])

    if debug:
        print('Sampled hyperparameters are {} {}'.format(rank_1, rank_2))
        print('Uniform map coordinates are {} {}'.format(coordinate_1, coordinate_2))
        print('Sampled realisation is {}'.format(index))
        print(uu_1)
        print(uu_2)

    # Write info to datablock
    block[pz, 'nbin'] = nbins
    block[pz, 'z'] = z
    block[pz, 'nz'] = len(z)

    for ibin in np.arange(nbins):
        block[pz, 'bin_{0}'.format(ibin+1)] = nz_sampled[ibin]
        block['ranks', 'mean_z_{0}'.format(ibin)] = nz_mean[index, ibin]

    return 0


def cleanup(config):

    return 0
