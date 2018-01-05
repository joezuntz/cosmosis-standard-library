from __future__ import print_function
from builtins import range
from builtins import object
import numpy as np
from numpy import pi
import os

# Some constants for this data set
D3000 = 3000 * 3001 / (2 * pi)
BETA = 0.0012309
DIPOLE_COSINE = -0.4033
NBIN = 56
NBAND = 2
BEAM_TERMS = 2


class SPTPolTheoryModel(object):
    def __init__(self, beamerr_file, windows_dir, use_te=True, use_ee=True):
        self.use_te = use_te
        self.use_ee = use_ee
        self.beam_errors = self.load_beam_errors(beamerr_file)
        self.ell, self.windows = self.load_windows(windows_dir)
        # self.r is the conversion from C_ell to D_ell
        self.r = self.ell * (self.ell + 1) / (2 * pi)

    def nuisance_parameter_names(self):
        p = ['t_cal', 'p_cal', 'kappa', 'beam_1', 'beam_2']
        if self.use_te:
            p += ['poisson_te', 'A_dust_te', 'alpha_dust_te']
        if self.use_ee:
            p += ['poisson_ee', 'A_dust_ee', 'alpha_dust_ee']
        return p

    def load_windows(self, windows_dir):
        # Load the first window just to check the ell values out
        filename = os.path.join(windows_dir, "window_{}".format(1))
        print("Loading band-power windows from {} directory".format(windows_dir))
        ell = np.loadtxt(filename)[:, 0].astype(int)
        n_ell = len(ell)
        n_band = 0
        self.band_indices = {}
        if self.use_te and self.use_ee:
            nband = 2
            ee_band = 0
            te_band = 1
            self.band_indices['te'] = te_band
            self.band_indices['ee'] = ee_band
        elif self.use_te:
            nband = 1
            te_band = 0
            self.band_indices['te'] = te_band
        elif self.use_ee:
            nband = 1
            ee_band = 0
            self.band_indices['ee'] = ee_band
        else:
            raise ValueError(
                "Must use at least one of te, ee in SPTPol likelihood")

        ee_start = 1
        te_start = ee_start + NBIN

        W = np.zeros((NBIN, n_ell, nband))

        if self.use_ee:
            for i in range(NBIN):
                filename = os.path.join(
                    windows_dir, "window_{}".format(ee_start + i))
                _, W_i = np.loadtxt(filename).T
                W[i, :, ee_band] = W_i

        if self.use_te:
            for i in range(NBIN):
                filename = os.path.join(
                    windows_dir, "window_{}".format(te_start + i))
                _, W_i = np.loadtxt(filename).T
                W[i, :, te_band] = W_i

        return ell, W

    def load_beam_errors(self, beamerr_file):
        print("Loading beam errors from file {}".format(beamerr_file))

        _, errors = np.loadtxt(beamerr_file).T
        beam_errors = {}
        if self.use_ee:
            beam_errors['ee'] = np.split(
                errors[:BEAM_TERMS * NBIN], BEAM_TERMS)
        if self.use_te:
            beam_errors['te'] = np.split(
                errors[BEAM_TERMS * NBIN:], BEAM_TERMS)
        return beam_errors

    def poisson(self, spectrum, cmb_theory_cl, nuisance_parameters):
        # The Poisson contribution ~ ell**2
        # Since the D_ell (which are what this function returns) already have the ell(ell+1)
        # quantity in them then we just use that conversion, except we normalize to ell at 3000.
        # i.e. we are approximating here that ell(ell+1)==ell**2
        # There is one parameter per spectrum (EE/TE) normalizing it.
        D = nuisance_parameters['poisson_' + spectrum] * self.r / D3000
        return D

    def supersample_lensing(self, spectrum, cmb_theory_cl, nuisance_parameters):
        # We have two different sets of ell here - the sample ell values (self.ell)
        # and all the points at which the theory is defined.  The latter must be a wider range
        # because we need the derivative of C_ell, which we do with finite differences
        ell = self.ell
        ell_all = cmb_theory_cl['ell']

        # ell**2 C_ell quantity
        C2 = cmb_theory_cl[spectrum] * ell_all**2

        # Deriviative d(ell**2 C_ell) / d(ell).
        # We are doing this with finite differences here - (Delta C_ell) / (Delta_ell) where
        # Delta_ell = 2
        derivative = 0.5 * (C2[ell + 1] - C2[ell - 1])

        # Total C_ell quantities from supersample lensing
        kappa = nuisance_parameters['kappa']
        C = -kappa * derivative / ell**2

        # Return them converted to D_ell = ell*(ell+1)*C_ell/2pi
        D = self.r * C
        return D

    def aberration(self, spectrum, cmb_theory_cl, nuisance_parameters):
        # ell sample values
        ell = self.ell

        # C_ell quantity
        C = cmb_theory_cl[spectrum]

        # C_ell derivative - again, this is a finite difference
        # but this time just on C_ell not ell**2 C_ell
        derivative = 0.5 * (C[ell + 1] - C[ell - 1])

        # Overall aberration C_ell contribution
        cl_ab = derivative * BETA * DIPOLE_COSINE / ell

        # return D_ell
        D = self.r * cl_ab
        return D

    def dust(self, spectrum, cmb_theory_cl, nuisance_parameters):
        # This is a two parameter model for each of spectrum=[TE,EE].
        # For the dust model D_ell ~ ell**2 (i.e. C_ell ~ ell**4)

        # The amplitude parameter A
        A = nuisance_parameters['A_dust_' + spectrum]

        # And the index parameter alpha
        alpha = nuisance_parameters['alpha_dust_' + spectrum]

        # The total D_ell spectrum
        D = A * (self.ell / 80.)**alpha
        return D

    def bin_into_bandpowers(self, spectrum, D):
        i = self.band_indices[spectrum]
        W = self.windows[:, :, i]
        B = np.einsum("bl,l->b", W, D)
        return B

    def apply_calibration(self, spectrum, B, nuisance_parameters):
        cal = nuisance_parameters['t_cal']**2
        if spectrum == "te":
            cal *= nuisance_parameters['p_cal']
        elif spectrum == "ee":
            cal *= nuisance_parameters['p_cal']**2
        return B / cal

    def apply_beam(self, spectrum, B, nuisance_parameters):
        # The beam described in the paper is different to the one
        # in the cosmomc code, though they are the same at first order
        # and the correction will be small.  Here we use the version
        # in the paper.
        beam = 1
        for i in range(BEAM_TERMS):
            f = nuisance_parameters["beam_{}".format(i + 1)]
            beam += f * self.beam_errors[spectrum][i]
        return B * beam

    def theory_plot(self, D0, D1, D2, D3, D4, backend="agg"):
        import matplotlib
        matplotlib.use(backend)
        import matplotlib.pyplot as plt
        plt.figure()
        plt.loglog(self.ell, abs(D0), label="CMB")
        plt.loglog(self.ell, abs(D1), label="Poisson")
        plt.loglog(self.ell, abs(D2), label="Supersample lensing")
        plt.loglog(self.ell, abs(D3), label="Aberration")
        plt.loglog(self.ell, abs(D4), label="Dust")
        plt.legend()
        plt.savefig(spectrum + ".png")
        plt.close()

    def bandpower_theory_model(self, spectrum, cmb_theory_cl, nuisance_parameters):
        # The various components
        D0 = cmb_theory_cl[spectrum][self.ell]
        D1 = self.poisson(spectrum, cmb_theory_cl, nuisance_parameters)
        D2 = self.supersample_lensing(
            spectrum, cmb_theory_cl, nuisance_parameters)
        D3 = self.aberration(spectrum, cmb_theory_cl, nuisance_parameters)
        D4 = self.dust(spectrum, cmb_theory_cl, nuisance_parameters)

        # The total predicted signal D_ell
        D = D0 + D1 + D2 + D3 + D4

        # We also return the components, in case interesting.
        components = {"ell": self.ell, "cmb": D0, "poisson": D1, "supersample_lensing": D2,
                      "aberration": D3, "dust": D4}

        # Converting to band powers B and applying multiplicative effects
        B = self.bin_into_bandpowers(spectrum, D)
        B = self.apply_calibration(spectrum, B, nuisance_parameters)
        B = self.apply_beam(spectrum, B, nuisance_parameters)

        # Return both
        return B, components

    def __call__(self, cmb_theory_cl, nuisance_parameters):
        """
        cmb_theory_cl: dictionary with keys ell, te, ee to numpy arrays of ell and spectra
        nuisance_parameters: dictionary of nuisance parameters
        """
        B = []

        if self.use_te:
            B_te, components_te = self.bandpower_theory_model(
                "te", cmb_theory_cl, nuisance_parameters)
            B.append(B_te)

        if self.use_ee:
            B_ee, components_ee = self.bandpower_theory_model(
                "ee", cmb_theory_cl, nuisance_parameters)
            B.append(B_ee)

        # Collect together complete data vector
        B = np.concatenate(B)
        return B, components_te, components_ee


class SPTPolData(object):
    def __init__(self, data_file, covmat_file, use_te=True, use_ee=True):
        self.use_te = use_te
        self.use_ee = use_ee
        self.vector = self.load_vector(data_file)
        self.covmat = self.load_covmat(covmat_file)

    def load_vector(self, data_file):
        print("Loading data vector from {}".format(data_file))
        data = np.loadtxt(data_file).T[1]
        vectors = []
        # For reasons that pass all understanding, the SPTpol data set supplied
        # a data vector that also includes the TT data, with the ordering
        # TE,EE,TT
        if self.use_te:
            vectors.append(data[:NBIN])
        if self.use_ee:
            vectors.append(data[NBIN:2 * NBIN])
        return np.concatenate(vectors)

    def load_covmat(self, covmat_file):
        print("Loading covariance matrix from {}".format(covmat_file))
        data = np.loadtxt(covmat_file).T.reshape((NBIN * 2, NBIN * 2))

        # Pull out the chunks we want.
        # The covariance matrix does *not* match the data vector supplied;
        # they have a different components - the covmat is just TE,EE (no TT).
        if self.use_te and self.use_ee:
            covmat = data
        elif self.use_te:
            covmat = data[:NBIN, :NBIN]
        elif self.use_ee:
            covmat = data[NBIN:, NBIN:]

        return covmat
