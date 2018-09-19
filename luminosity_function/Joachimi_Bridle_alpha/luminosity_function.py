from builtins import range
import numpy as np
from scipy import interpolate


def jb_calculate_alpha(a, zmax, Nz):
    # First set up a redshift array for the given survey parameters
    z = np.linspace(0, zmax, Nz)

    # Then just apply the coefficients
    alpha = np.zeros_like(z)
    alpha += a[0]
    alpha += a[1] * z
    alpha += a[2] * z**2
    return alpha, z


def initialise_jb_coefficients(mag_lim):
    b = np.array([[0.44827, 	0.0, 		0.0],
                  [-1.0, 		1.0, 		1.0],
                  [0.05617, 	0.19658,	0.18107],
                  [0.07704, 	3.31359,	3.05213],
                  [-11.3768,	-2.5028,	-2.5027]])

    a = np.array([b[0][0] + b[1][0] * (b[2][0] * mag_lim - b[3][0])**b[4][0],
                  b[0][1] + b[1][1] * (b[2][1] * mag_lim - b[3][1])**b[4][1],
                  b[0][2] + b[1][2] * (b[2][2] * mag_lim - b[3][2])**b[4][2]])

    return a


def get_binned_alpha(block, alpha, z, sample='wl_number_density'):
    n_z, z1 = load_n_z(block, sample)

    z_med = evaluate_median_z(n_z, z1)

    interpolator = interpolate.interp1d(z, alpha)
    alpha_binned = interpolator(z_med)

    return alpha_binned, z_med


def load_n_z(block, sample='wl_number_density'):
    """ Load the n(z) profile in each bin as a 2d array """
    num_den = sample

    n_z = []
    N_zbins = block[num_den, 'nbin']
    for i in range(1, N_zbins + 1):
        n_z += [block.get_double_array_1d(num_den, 'bin_%d' % i)]
    n_z = np.array(n_z)

    z = block.get_double_array_1d(num_den, 'z')

    return n_z, z


def evaluate_mean_z(n_z, z):
    """Integrate n(z) in each bin to get the mean redshift. """
    z_mean = []
    for n in n_z:
        z_mean += [sum(n * z) / sum(n)]

    z_mean = np.array(z_mean)

    return z_mean


def evaluate_median_z(n_z, z):
    """Evaluate the median of the redshift distribution for each bin. """
    import scipy.stats as stats
    z_med = []
    for n in n_z:
        prob_dist = n / n.sum()
        gen = stats.rv_discrete(values=(z, prob_dist), inc=z[1] - z[0])
        z_med += [gen.median()]
    z_med = np.array(z_med)

    return z_med
