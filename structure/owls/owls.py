"""
This module loads data from the powtable files which summarize the OWLS
results made by Tim Eifler et al.

I interpolate into that data using a bivariate spline to get an estimate of the
effect on the matter power from baryons at a given z and k.

This requires an additional parameter in the range (0,1) which controls how far between the 
two extremal feedback scenarios we live.

I could have made mistakes here easily in understanding what is happening! In particular
I'm not totally sure the value of the read-in data should be squared.

"""

from builtins import range
from builtins import object
import numpy as np
try:
    import scipy.interpolate
except ImportError:
    sys.stderr.write(
        "The OWLS baryon power code requires the scipy python module\n")
    sys.stderr.write("but it could not be found.  This code will fail.\n")
    raise ImportError(
        "Require scipy module to run baryon code.  Could not find it.")

import os
dirname = os.path.split(__file__)[0]

# File names
DM_FILENAME = os.path.join(dirname, 'powtable_DMONLY_all.dat')
UPPER_FILENAME = os.path.join(dirname, 'powtable_NOSN_all.dat')
LOWER_FILENAME = os.path.join(dirname, 'powtable_AGN_all.dat')

MIN_K_MODULATION = 0.05
MIN_Z_MODULATION = 0.0
MAX_Z_MODULATION = 5.0


class OwlsFileUser(object):
    @staticmethod
    def _load_file(filename):
        data = np.loadtxt(filename, comments='%').T
        z = data[0]
        k = data[1]
        P = data[2]
        out_z = np.unique(data[0])
        out_z = out_z[out_z < 3.1]
        kmin = k.min()
        kmax = k.max()
        nk_out = 100
        nz_out = out_z.size
        out_k = np.linspace(np.log10(kmin), np.log10(kmax), nk_out)
        out_P = np.zeros((nk_out, nz_out))

        for i, z_i in enumerate(out_z):
            w = (z == z_i)
            k_i = k[w]
            P_i = P[w]
            # if z_i==0.5:
            # 	import pylab
            # 	pylab.loglog(k_i, P_i, label=filename)
            ok = P_i > 0
            k_i = np.log10(k_i[ok])
            P_i = np.log10(P_i[ok])
            S = scipy.interpolate.InterpolatedUnivariateSpline(
                k_i, P_i)  # , s=3)
            P_fit = S(out_k)
            out_P[:, i] = P_fit

        return out_k, out_z, out_P


class BaryonPowerModulator(OwlsFileUser):
    """Object containing stored data from baryon power, with a method to modulate an input P(k,z) grid according to a single parameter"""

    def __init__(self, dark_matter_filename=DM_FILENAME,
                 upper_filename=UPPER_FILENAME,
                 lower_filename=LOWER_FILENAME):
        # Load the DM-only file, which forms the baseline.
        # Generate an interpolator from it.
        super(BaryonPowerModulator, self).__init__()
        dm_k, dm_z, dm_P = self._load_file(dark_matter_filename)
        self.dm = scipy.interpolate.RectBivariateSpline(dm_k, dm_z, dm_P)
        norm = dm_P[0, 0]  # normalize at z=0, k=smallest

        # Now do the same for the upper limit, but divide by the value
        # for DM-only to get the ratio
        u_k, u_z, u_P = self._load_file(upper_filename)
        u_P *= norm / u_P[0, 0]
        self.upper = scipy.interpolate.RectBivariateSpline(u_k, u_z, u_P)

        l_k, l_z, l_P = self._load_file(lower_filename)
        l_P *= norm / l_P[0, 0]
        self.lower = scipy.interpolate.RectBivariateSpline(l_k, l_z, l_P)

        # import pylab
        # z = 1.0
        # k = np.log10(np.logspace(-1,3, 100))
        # L = 10**self.lower(k,z)
        # U = 10**self.upper(k,z)
        # M = 10**self.dm(k, z)
        # pylab.loglog(10**k, L)
        # pylab.loglog(10**k, U)
        # pylab.loglog(10**k, M)
        # pylab.semilogx(10**k, (L-U)/L,','	)
        # pylab.legend()
        # pylab.show()
        # jkhgjhghj

    def modulate(self, k, z, P, r):
        logk = np.log10(k)
        d = self.dm(logk, z)
        u = self.upper(logk, z) - d
        l = self.lower(logk, z) - d
        modulation = r * 10**u.squeeze() + (1 - r) * 10**l.squeeze()
        return P * modulation


class ChebyshevBaryonPowerModulator(BaryonPowerModulator):
    def __init__(self, *args, **kwargs):
        self.extremes = kwargs.pop("extremes", 0.0)
        self.nterm = kwargs.pop("nterm")
        super(ChebyshevBaryonPowerModulator, self).__init__(*args, **kwargs)

    @staticmethod
    def _chebyshev(x, n):
        return np.cos(n * np.arccos(x))

    def modulate(self, k, z, P, r):
        """ r_z and r_logk should be in 0,1"""
        r_k = r[:self.nterm]
        r_z = r[self.nterm:]
        assert len(r_k) == len(
            r_z) == self.nterm, "Parameters not found correctly for OWLS code"
        logk = np.log10(k)
        d = self.dm(logk, z)
        u = 10**(self.upper(logk, z) - d)
        l = 10**(self.lower(logk, z) - d)
        u = u.squeeze()
        l = l.squeeze()
        # we now have (l,u) into which we want to interpolate
        kmin = np.log10(MIN_K_MODULATION)
        kmax = logk.max()
        zmin = MIN_Z_MODULATION
        zmax = MAX_Z_MODULATION

        n_zcoeffs = len(r_z)
        n_kcoeffs = len(r_k)
        # logk and z, scaled to (0,1)
        ks = (logk - kmin) / (kmax - kmin)
        zs = (z - zmin) / (zmax - zmin)

        # move the coefficients to range (-1,+1)
        k_coeffs = 2 * np.array(r_k) - 1
        z_coeffs = 2 * np.array(r_z) - 1

        # get the k and z based
        Tk = [self._chebyshev(ks, i) for i in range(n_kcoeffs)]
        Tz = [self._chebyshev(zs, i) for i in range(n_zcoeffs)]
        Tk = np.array(Tk)
        Tz = np.array(Tz)
        pk = np.dot(k_coeffs, Tk)
        pz = np.dot(z_coeffs, Tz)

        # Ignore things outside our range
        pk[logk < kmin] = 0.0
        pz[z > zmax] = 0.0

        # pk and pz are both in range (-1,1)
        # so p will be too.
        # bring it to the range (0,1)
        # with some leeway given by the extremeness parameters
        p = np.outer(pk, pz)
        p = (p + 1) / 2
        p = p.clip(-self.extremes, 1 + self.extremes)

        modulation = p * u + (1 - p) * l
        return P * modulation


class FixedBaryonPowerModulator(OwlsFileUser):
    """ Single value baryon marginalization based on a loaded OWLS table"""

    def __init__(self, simulation_filename, dark_matter_filename=DM_FILENAME):
        super(FixedBaryonPowerModulator, self).__init__()
        dm_k, dm_z, dm_P = self._load_file(dark_matter_filename)
        self.dm = scipy.interpolate.RectBivariateSpline(dm_k, dm_z, dm_P)
        norm = dm_P[0, 0]  # normalize at z=0, k=smallest

        sim_k, sim_z, sim_P = self._load_file(simulation_filename)
        self.sim_k = sim_k
        # sim_P[0,0]=dm_P[0,0]
        #sim_P *= norm/sim_P[0,0]
        self.sim = scipy.interpolate.RectBivariateSpline(sim_k, sim_z, sim_P)

    def modulate(self, k, z, P, r=None):
        # The r argument means we have the same calling structure as the other ones
        # It is ignored
        logk = np.log10(k)
        s = (10**self.sim(logk, z) - 10**self.dm(logk, z)) / \
            (10**self.dm(logk, z))
        s[k < self.sim_k[0]] = 0.
        modulation = 1 + s.squeeze()
        return P * modulation


class ScaledBaryonPowerModulator(OwlsFileUser):
    """ Scaled baryon modulation based on a loaded OWLS table"""

    def __init__(self, simulation_filename, dark_matter_filename=DM_FILENAME):
        super(ScaledBaryonPowerModulator, self).__init__()
        dm_k, dm_z, dm_P = self._load_file(dark_matter_filename)
        self.dm = scipy.interpolate.RectBivariateSpline(dm_k, dm_z, dm_P)
        norm = dm_P[0, 0]  # normalize at z=0, k=smallest

        sim_k, sim_z, sim_P = self._load_file(simulation_filename)
        self.sim_k = sim_k
        sim_P[0, 0] = dm_P[0, 0]
        #sim_P *= norm/sim_P[0,0]
        self.sim = scipy.interpolate.RectBivariateSpline(sim_k, sim_z, sim_P)

    def modulate(self, k, z, P, r=None):
        logk = np.log10(k)
        s = (10**self.sim(logk, z) - 10**self.dm(logk, z)) / \
            (10**self.dm(logk, z))
        s[k < self.sim_k[0]] = 0.
        modulation = 1 + r * s.squeeze()
        return P * modulation
