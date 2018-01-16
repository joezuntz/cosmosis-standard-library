from builtins import range
import numpy as np
from scipy.interpolate import interp1d
import warnings


def compute_gaussian_covariance(sky_area, get_theory_cl, block, AB, CD):
    """ 
    Calculates an analytic Gaussian covariance matrix for a  
    given combination of spectra and redshift bins using the 
    theory Cls. See e.g. Joachimi and Bridle (2010) eq. (35),
    whose notation we follow.

    The covariance is between C_ell^{AB} and C_ell^{CD} for bin pairs (i,j)
    and (k,l) respectively.

    AB and CD must be shear-shear, shear-position, or position-position;
    we have not yet coded the other terms.

    In this code AB and CD are SpectrumMeasurement objects. ij and kl
    are pairs of integers.
    """

    # These are all one of
    # galaxy_shear_emode_fourier, galaxy_position_fourier
    A = AB.type1.name
    B = AB.type2.name
    C = CD.type1.name
    D = CD.type2.name

    # Check that these are the kind of spectra this code is ready for.
    for X in (A, B, C, D):
        assert X in ["galaxy_shear_emode_fourier",
                     "galaxy_position_fourier"], "Have not yet coded up gaussian covariance for spectra other than pos-pos, shear-shear, shear-pos"

    # Data in the SpectrumMeasurement is stored as a set of vectors specifying
    # the two tomorgraphic bins each data point refers to, and also the angular
    # value ell.  Each of these can be different for the two spectra. For example
    # we might not have cross-spectra for position-position but we probably will
    # for lensing.
    ell_AB = AB.angle
    ell_CD = CD.angle

    # These are the tomographic bin indices
    bin1_AB = AB.bin1
    bin2_AB = AB.bin2
    bin1_CD = CD.bin1
    bin2_CD = CD.bin2

    # The total number of data points for each of the spectra.
    n_AB = len(ell_AB)
    n_CD = len(ell_CD)

    # This will be our actual output matrix
    covmat = np.zeros((n_AB, n_CD))

    # There might be no ell values in one of the
    # bins if one has been masked out or something.
    # This is fine; we will handle it in the calling function.
    if n_AB == 0 or n_CD == 0:
        return covmat

    # Work out the delta ells for each measurement.
    delta_ell_AB = compute_delta_ells(AB)
    delta_ell_CD = compute_delta_ells(CD)

    # Now loop through all the places where we are non-zero.
    # The find_equal_ell thing checks if they are the same and returns
    # the indices.
    for x, y, ell in find_equal_ell(ell_AB, ell_CD):

        i = bin1_AB[x]
        j = bin2_AB[x]
        k = bin1_CD[y]
        l = bin2_CD[y]

        # Now we can look the various C_ell values we want using our
        # lookup function. This is where things will fail if, for example,
        # we did not predict NE but want the covariance between NN and EE.
        # This function (which is passed in as an argument) could just do
        # a simple lookup in the block but might also use cached values
        # as in the case we have in 2pt_like where the values have already
        # been loaded
        C_AC_ik = get_theory_cl(block, A, C, i, k, ell)
        C_BD_jl = get_theory_cl(block, B, D, j, l, ell)
        C_AD_il = get_theory_cl(block, A, D, i, l, ell)
        C_BC_jk = get_theory_cl(block, B, C, j, k, ell)

        # This is a reasonable hack. Get the delta_ell for each of the values
        # independently and take their geometric mean. They should really
        # just be equal for all this to make any sense.
        #delta_ell = np.sqrt(delta_ell_AB[(i,j)](ell)*delta_ell_CD[(i,j)](ell))
        # If this bin is not here just get any of them.
        delta_AB = delta_ell_AB.get((i, j))
        delta_CD = delta_ell_CD.get((i, j))
        if delta_AB is None:
            warnings.warn(
                "There are no delta-ell values for some combinations of ({},{}).".format(A, B))
            delta_AB = list(delta_ell_AB.values())[0]
        if delta_CD is None:
            warnings.warn(
                "There are no delta-ell values for some combinations ({},{}).".format(C, D))
            delta_CD = list(delta_ell_CD.values())[0]

        delta_AB = delta_AB(ell)
        delta_CD = delta_CD(ell)

        delta_ell = np.sqrt(delta_AB * delta_CD)

        # The pre-factor.
        p = 2. * np.pi / (ell * delta_ell * sky_area)
        covmat[x, y] = p * (C_AC_ik * C_BD_jl + C_AD_il * C_BC_jk)

    return covmat


def find_equal_ell(ell_1, ell_2):
    n_1 = len(ell_1)
    n_2 = len(ell_2)
    for i in range(n_1):
        ell_i = ell_1[i]
        for j in range(n_2):
            if ell_i == ell_2[j]:
                yield i, j, ell_i

# thanks to:
# http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range


def extrap1d(interpolator):
    "Turn an interpolator into a (linear) extrapolator"
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        else:
            return interpolator(x)

    return pointwise


def compute_delta_ells(XY):
    """This little utility function is used in the gaussian covariance code
    to collect delta-ell values from a spectrum
    """
    output = {}
    for i in set(XY.bin1):
        for j in set(XY.bin2):
            ell, _ = XY.get_pair(i, j)
            if len(ell) == 0:
                continue

            d = np.diff(ell)
            # add another copy of the first diff value
            delta = np.concatenate([d[:1], d])
            output[(i, j)] = extrap1d(interp1d(ell, delta))
    return output
