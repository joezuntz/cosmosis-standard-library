from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
import scipy.interpolate
import pdb
import scipy.stats.mstats as mstats
import numpy as np


class Cl_class(object):
    def __init__(self, block, config):

        # Spectra to use
        self.shear = config['shear']
        self.intrinsic_alignments = (
            config['intrinsic_alignments'], config['GI'], config['II'])
        self.position = config['position']
        self.ggl = config['ggl']
        self.magnification = config['magnification']
        self.cmb_kappa = config['cmb_kappa']
        self.kappa_shear = config['kappa_shear']
        self.kappa_pos = config['kappa_pos']

        self.noise = config['noise']
        self.bias = config['bias'][0]
        self.m_per_bin = config['bias'][1]

        shear_cat = config['shear_cat']
        pos_cat = config['pos_cat']

        self.dobinning = config['binning']
        self.window = config['window']

        # Relevant parameters for the noise
        if self.noise:
            self.sigma_gamma = config['shape_dispersion']
            self.ngal_shear, self.ngal_pos = config['ngal']

        # And the configuration of the theory spectra
        self.Nzbin_shear = int(block.get_int(shear_cat, 'nzbin', default=0))
        self.Nzbin_pos = int(block.get_int(shear_cat, 'nzbin', default=0))
        self.Nzbin_cmb = 1.
        if self.shear:
            self.zbin_edges_shear = [block[shear_cat, 'edge_%d' % i]
                                     for i in range(1, self.Nzbin_shear + 1)]
            if self.intrinsic_alignments[0]:
                gg_section = 'shear_cl_gg'
            else:
                gg_section = 'shear_cl'
            self.l_shear = block[gg_section, 'ell']
            self.Nl_shear = len(self.l_shear)
            if self.dobinning:
                self.Nlbin_shear = int(config['nlbin_shear'])
            else:
                self.Nlbin_shear = self.Nl_shear
        if self.position:
            self.zbin_edges_pos = [block[pos_cat, 'edge_%d' % i]
                                   for i in range(1, self.Nzbin_pos + 1)]
            self.l_pos = block['matter_cl', 'ell']
            self.Nl_pos = len(self.l_pos)
            if self.dobinning:
                self.Nlbin_pos = int(config['nlbin_pos'])
            else:
                self.Nlbin_pos = self.Nl_pos
        if self.ggl:
            self.l_ggl = block['matter_cl', 'ell']
            self.Nl_ggl = len(self.l_ggl)
            if self.dobinning:
                self.Nlbin_ggl = int(config['nlbin_ggl'])
            else:
                self.Nlbin_ggl = self.Nl_ggl

        if self.cmb_kappa:
            self.l_kk = block['cmb_kappa_cl', 'ell']
            self.Nl_kk = len(self.l_pos)
            if self.dobinning:
                self.Nlbin_kk = int(config['nlbin_kk'])
            else:
                self.Nlbin_kk = self.Nl_kk
        if self.kappa_shear:
            self.l_ke = block['cmb_kappa_shear_cl', 'ell']
            self.Nl_ke = len(self.l_pos)
            if self.dobinning:
                self.Nlbin_ke = int(config['nlbin_ke'])
            else:
                self.Nlbin_ke = self.Nl_ke

        if self.kappa_pos:
            self.l_kn = block['cmb_kappa_matter_cl', 'ell']
            self.Nl_kn = len(self.l_kn)
            if self.dobinning:
                self.Nlbin_kn = int(config['nlbin_kn'])
            else:
                self.Nlbin_kn = self.Nl_kn

        if self.bias:
            self.multiplicative_bias = [block[shear_cat, "m%d" % i]
                                        for i in range(1, self.Nzbin_shear + 1)]

        # Finally get the desired binning
        self.get_l_bins(config)

    def load_and_generate_observable_cls(self, block, names):

        # Set up somewhere to put the observable spectra
        if self.shear:
            self.C_ee = np.zeros(
                (self.Nzbin_shear, self.Nzbin_shear, self.Nl_shear))
            self.C_ee_binned = np.zeros(
                (self.Nzbin_shear, self.Nzbin_shear, self.Nlbin_shear))
        if self.position:
            self.C_nn = np.zeros((self.Nzbin_pos, self.Nzbin_pos, self.Nl_pos))
            self.C_nn_binned = np.zeros(
                (self.Nzbin_pos, self.Nzbin_pos, self.Nlbin_pos))
        if self.ggl:
            self.C_ne = np.zeros(
                (self.Nzbin_pos, self.Nzbin_shear, self.Nl_pos))
            self.C_ne_binned = np.zeros(
                (self.Nzbin_pos, self.Nzbin_shear, self.Nlbin_ggl))

        # Then cycle through all the redshift bin combinations
        if self.shear:
            for i in range(1, self.Nzbin_shear + 1):
                for j in range(1, self.Nzbin_shear + 1):
                    bin = "bin_%d_%d" % (i, j)
                    bin_tr = "bin_%d_%d" % (j, i)

                    # The C_GG,II,mm,gg spectra are symmetric
                    # This is just bookkeeping to account for the fact we only have half of them
                    if (j < i):
                        a = bin
                    else:
                        a = bin_tr
                        # GG
                    if self.intrinsic_alignments[0]:
                        self.C_ee[i - 1][j -
                                         1] += block.get_double_array_1d("shear_cl_gg", a)
                        if self.intrinsic_alignments[1]:
                            # GI
                            self.C_ee[i - 1][j -
                                             1] += block.get_double_array_1d("shear_cl_gi", bin)
                            # IG
                            self.C_ee[i - 1][j -
                                             1] += block.get_double_array_1d("shear_cl_gi", bin_tr)
                        if self.intrinsic_alignments[2]:
                            # II
                            self.C_ee[i - 1][j -
                                             1] += block.get_double_array_1d("shear_cl_ii", a)
                    else:
                        self.C_ee[i - 1][j -
                                         1] += block.get_double_array_1d("shear_cl", a)
            block["galaxy_shape_cl_unbinned", "ell"] = block.get_double_array_1d(
                "shear_cl_gg", "ell")

        if self.position:
            for i in range(1, self.Nzbin_pos + 1):
                for j in range(1, self.Nzbin_pos + 1):
                    bin = "bin_%d_%d" % (i, j)
                    bin_tr = "bin_%d_%d" % (j, i)

                    # The C_GG,II,mm,gg spectra are symmetric
                    # This is just bookkeeping to account for the fact we only have half of them
                    if (j < i):
                        a = bin
                    else:
                        a = bin_tr

                    # gg
                    self.C_nn[i - 1][j -
                                     1] += block.get_double_array_1d('matter_cl', a)
                    if self.magnification:
                        # mg
                        self.C_nn[i - 1][j - 1] += block.get_double_array_1d(
                            "galaxy_magnification_cl", bin)
                        self.C_nn[i - 1][j - 1] += block.get_double_array_1d(
                            "galaxy_magnification_cl", bin_tr)  # gm
                        self.C_nn[i - 1][j - 1] += block.get_double_array_1d(
                            "magnification_magnification_cl", a)  # mm

            block["galaxy_position_cl_unbinned",
                  "ell"] = block.get_double_array_1d("matter_cl", "ell")
        if self.ggl:
            block["galaxy_position_shape_cross_cl_unbinned",
                  "ell"] = block.get_double_array_1d("matter_cl", "ell")
        for i in range(1, self.Nzbin_pos + 1):
            for j in range(1, self.Nzbin_shear + 1):
                bin = "bin_%d_%d" % (i, j)
                bin_tr = "bin_%d_%d" % (j, i)

                # The C_GG,II,mm,gg spectra are symmetric
                # This is just bookkeeping to account for the fact we only have half of them
                if (j < i):
                    a = bin
                else:
                    a = bin_tr

                if self.ggl:
                    # gG
                    self.C_ne[i - 1][j -
                                     1] += block.get_double_array_1d("ggl_cl", bin)
                    if self.intrinsic_alignments[0]:
                        # gI
                        self.C_ne[i - 1][j -
                                         1] += block.get_double_array_1d("gal_IA_cross_cl", bin)
                        if self.magnification:
                            self.C_ne[i - 1][j - 1] += block.get_double_array_1d(
                                "magnification_intrinsic_cl", bin)  # mI
                    if self.magnification:
                        # mG
                        self.C_ne[i - 1][j - 1] += block.get_double_array_1d(
                            "magnification_shear_cl", bin)

                if not self.noise:
                    # Finally resample the spectra in the survey angular frequency bins
                    if self.shear:
                        self.C_ee_binned[i - 1][j - 1] = get_binned_cl(self.C_ee[i - 1][j - 1], self.l_shear,
                                                                       self.lbin_edges_shear, self.dobinning, self.window)
                        if self.bias:
                            self.apply_measurement_bias(i, j, 'shear')
                        block["galaxy_shape_cl_unbinned",
                              a] = self.C_ee[i - 1][j - 1]
                    if self.position:
                        self.C_nn_binned[i - 1][j - 1] = get_binned_cl(self.C_nn[i - 1][j - 1], self.l_pos,
                                                                       self.lbin_edges_pos, self.dobinning, self.window)
                        block["galaxy_position_cl_UNBINNED",
                              a] = self.C_nn[i - 1][j - 1]
                    if self.ggl:
                        self.C_ne_binned[i - 1][j - 1] = get_binned_cl(self.C_ne[i - 1][j - 1], self.l_ggl,
                                                                       self.lbin_edges_pos, self.dobinning, self.window)
                        if self.bias:
                            self.apply_measurement_bias(i, j, 'ggl')
                        block["galaxy_position_shape_cross_cl_unbinned",
                              a] = self.C_ne[i - 1][j - 1]
        if self.noise:
            # Add shot noise if required
            self.add_noise(block)
            # If noise was added earlier, the binning is done here rather than
            # immediately on loading
            if self.shear:
                for i in range(1, self.Nzbin_shear + 1):
                    for j in range(1, self.Nzbin_shear + 1):
                        self.C_ee_binned[i - 1][j - 1] = get_binned_cl(self.C_ee[i - 1][j - 1], self.l_shear,
                                                                       self.lbin_edges_shear, self.dobinning, self.window)
                        if self.bias:
                            self.apply_measurement_bias(i, j, 'shear')
                        block["galaxy_shape_cl_unbinned", "bin_%d_%d" %
                              (i, j)] = self.C_ee[i - 1][j - 1]
            if self.position:
                for i in range(1, self.Nzbin_pos + 1):
                    for j in range(1, self.Nzbin_pos + 1):
                        self.C_nn_binned[i - 1][j - 1] = get_binned_cl(self.C_nn[i - 1][j - 1], self.l_pos,
                                                                       self.lbin_edges_pos, self.dobinning, self.window)
                        block["galaxy_position_cl_UNBINNED", "bin_%d_%d" %
                              (i, j)] = self.C_nn[i - 1][j - 1]
                        print(a)
            if self.ggl:
                for i in range(1, self.Nzbin_pos + 1):
                    for j in range(1, self.Nzbin_shear + 1):
                        self.C_ne_binned[i - 1][j - 1] = get_binned_cl(self.C_ne[i - 1][j - 1], self.l_ggl,
                                                                       self.lbin_edges_pos, self.dobinning, self.window)
                        if self.bias:
                            self.apply_measurement_bias(i, j, 'ggl')
                        block["galaxy_position_shape_cross_cl_unbinned",
                              "bin_%d_%d" % (i, j)] = self.C_ne[i - 1][j - 1]

    def apply_measurement_bias(self, i, j, mode=None):
        if not self.m_per_bin:
            m0 = self.multiplicative_bias[0]

        # Compute scaling parameter for this pair of redshift bins
        if self.m_per_bin:
            mi = self.multiplicative_bias[i - 1]
            mj = self.multiplicative_bias[j - 1]
        else:
            mi, mj = m0, m0

        # Apply scaling
        if mode == 'shear':
            self.C_ee_binned[i - 1][j - 1] *= (1 + mi) * (1 + mj)
        if mode == 'ggl':
            self.C_ne_binned[i - 1][j - 1] *= (1 + mj)

    def add_noise(self, block):

        n_binned_shear = get_binned_number_densities(
            self.Nzbin_shear, self.ngal_shear)
        n_binned_pos = get_binned_number_densities(
            self.Nzbin_pos, self.ngal_pos)

        # Create noise matrices with the same shape as the Cls
        # These are diagonal in the x,z plane (fixed l) and constant along the y axis (constant redshift)

        N_ee_0 = np.identity(self.Nzbin_shear) * \
            self.sigma_gamma**2 / (2. * n_binned_shear)
        N_nn_0 = np.identity(self.Nzbin_pos) * 1. / n_binned_pos

        N_shot_ee = []
        N_shot_nn = []

        if self.shear:
            for i in range(len(self.C_ee[0][0])):
                N_shot_ee += [N_ee_0]
            N_shot_ee = np.swapaxes(N_shot_ee, 0, 2)
            N_shot_ee = np.swapaxes(N_shot_ee, 0, 1)
        if self.position:
            for i in range(len(self.C_nn[0][0])):
                N_shot_nn += [N_nn_0]

            N_shot_nn = np.swapaxes(N_shot_nn, 0, 2)
            N_shot_nn = np.swapaxes(N_shot_nn, 0, 1)

        # Then add the relevant noise to the Cl matrices
        if self.shear:
            self.C_ee += N_shot_ee
        if self.position:
            self.C_nn += N_shot_nn

    def get_l_bins(self, config):

        if self.dobinning:
            # Define some l bins for these galaxy samples
            lmin, lmax = config['lmin_shear'], config['lmax_shear']
            if self.shear:
                self.lbin_edges_shear = np.logspace(
                    np.log10(lmin), np.log10(lmax), self.Nlbin_shear + 1)
                self.l_bins_shear = np.zeros_like(self.lbin_edges_shear[:-1])
                if self.window == 'tophat':
                    self.l_bins_shear = np.exp(
                        (np.log(self.lbin_edges_shear[1:] * self.lbin_edges_shear[:-1])) / 2.0)
                elif self.window == 'tophat-arithmetic':
                    self.l_bins_shear = (
                        self.lbin_edges_shear[1:] + self.lbin_edges_shear[:-1]) / 2.0
                elif self.window == 'delta':
                    # Just take the mid point of each bin and sample the Cls there
                    for i in range(len(self.lbin_edges_shear) - 1):
                        lmin0 = self.lbin_edges_shear[i]
                        lmax0 = self.lbin_edges_shear[i + 1]
                        sel = (self.l_shear > lmin0) & (self.l_shear < lmax0)
                        l_in_window = self.l_shear[sel]
                        self.l_bins_shear[i] = l_in_window[len(
                            l_in_window) / 2]

            if self.position:
                lmin, lmax = config['lmin_pos'], config['lmax_pos']
                self.lbin_edges_pos = np.logspace(
                    np.log10(lmin), np.log10(lmax), self.Nlbin_pos + 1)
                self.l_bins_pos = np.zeros_like(self.lbin_edges_pos[:-1])
                if self.window == 'tophat':
                    self.l_bins_pos = np.exp(
                        (np.log(self.lbin_edges_pos[1:] * self.lbin_edges_pos[:-1])) / 2.0)
                elif self.window == 'tophat-arithmetic':
                    self.l_bins_pos = (
                        self.lbin_edges_pos[1:] + self.lbin_edges_pos[:-1]) / 2.0
                elif self.window == 'delta':
                    for i in range(len(self.lbin_edges_pos) - 1):
                        lmin0 = self.lbin_edges_pos[i]
                        lmax0 = self.lbin_edges_pos[i + 1]
                        sel = (self.l_pos > lmin0) & (self.l_pos < lmax0)
                        l_in_window = self.l_pos[sel]
                        self.l_bins_pos[i] = l_in_window[len(l_in_window) / 2]

            if self.position and self.shear:
                lmin, lmax = config['lmin_ggl'], config['lmax_ggl']
                self.lbin_edges_ggl = np.logspace(
                    np.log10(lmin), np.log10(lmax), self.Nlbin_ggl + 1)
                self.l_bins_ggl = np.zeros_like(self.lbin_edges_ggl[:-1])
                lmin, lmax = config['lmin_ggl'], config['lmax_ggl']
                self.lbin_edges_ggl = np.logspace(
                    np.log10(lmin), np.log10(lmax), self.Nlbin_ggl + 1)
                if self.window == 'tophat':
                    self.l_bins_ggl = np.exp(
                        (np.log(self.lbin_edges_ggl[1:] * self.lbin_edges_ggl[:-1])) / 2.0)
                elif self.window == 'tophat-arithmetic':
                    self.l_bins_ggl = (
                        self.lbin_edges_ggl[1:] + self.lbin_edges_ggl[:-1]) / 2.0
                elif self.window == 'delta':
                    for i in range(len(self.lbin_edges_ggl) - 1):
                        lmin0 = self.lbin_edges_ggl[i]
                        lmax0 = self.lbin_edges_ggl[i + 1]
                        sel = (self.l_ggl > lmin0) & (self.l_ggl < lmax0)
                        l_in_window = self.l_ggl[sel]
                        self.l_bins_ggl[i] = l_in_window[len(l_in_window) / 2]
            if self.cmb_kappa:
                lmin, lmax = config['lmin_cmb_kappa'], config['lmax_cmb_kappa']
                self.lbin_edges_kk = np.linspace(lmin, lmax, self.Nlbin_kk + 1)
                self.l_bins_kk = np.zeros_like(self.lbin_edges_kk[:-1])
                if self.window == 'tophat':
                    self.l_bins_kk = np.exp(
                        (np.log(self.lbin_edges_kk[1:] * self.lbin_edges_kk[:-1])) / 2.0)
                elif self.window == 'tophat-arithmetic':
                    self.l_bins_kk = (
                        self.lbin_edges_kk[1:] + self.lbin_edges_kk[:-1]) / 2.0
                elif self.window == 'delta':
                    for i in range(len(self.lbin_edges_kk) - 1):
                        lmin0 = self.lbin_edges_kk[i]
                        lmax0 = self.lbin_edges_kk[i + 1]
                        sel = (self.l_kk > lmin0) & (self.l_kk < lmax0)
                        l_in_window = self.l_kk[sel]
                        self.l_bins_kk[i] = l_in_window[len(l_in_window) / 2]

            if self.kappa_shear:
                lmin, lmax = config['lmin_kappa_shear'], config['lmax_kappa_shear']
                self.lbin_edges_ke = np.linspace(lmin, lmax, self.Nlbin_ke + 1)
                self.l_bins_ke = np.zeros_like(self.lbin_edges_ke[:-1])
                if self.window == 'tophat':
                    self.l_bins_ke = np.exp(
                        (np.log(self.lbin_edges_ke[1:] * self.lbin_edges_ke[:-1])) / 2.0)
                elif self.window == 'tophat-arithmetic':
                    self.l_bins_ke = (
                        self.lbin_edges_ke[1:] + self.lbin_edges_ke[:-1]) / 2.0
                elif self.window == 'delta':
                    for i in range(len(self.lbin_edges_ke) - 1):
                        lmin0 = self.lbin_edges_ke[i]
                        lmax0 = self.lbin_edges_ke[i + 1]
                        sel = (self.l_ke > lmin0) & (self.l_ke < lmax0)
                        l_in_window = self.l_ke[sel]
                        self.l_bins_ke[i] = l_in_window[len(l_in_window) / 2]
            if self.kappa_pos:
                lmin, lmax = config['lmin_kappa_pos'], config['lmax_kappa_pos']
                self.lbin_edges_kn = np.linspace(lmin, lmax, self.Nlbin_kn + 1)
                self.l_bins_kn = np.zeros_like(self.lbin_edges_kn[:-1])
                if self.window == 'tophat':
                    self.l_bins_kn = np.exp(
                        (np.log(self.lbin_edges_kn[1:] * self.lbin_edges_kn[:-1])) / 2.0)
                elif self.window == 'tophat-arithmetic':
                    self.l_bins_kn = (
                        self.lbin_edges_kn[1:] + self.lbin_edges_kn[:-1]) / 2.0
                elif self.window == 'delta':
                    for i in range(len(self.lbin_edges_kn) - 1):
                        lmin0 = self.lbin_edges_kn[i]
                        lmax0 = self.lbin_edges_kn[i + 1]
                        sel = (self.l_kn > lmin0) & (self.l_kn < lmax0)
                        l_in_window = self.l_kn[sel]
                        self.l_bins_kn[i] = l_in_window[len(l_in_window) / 2]
        else:
            if self.shear:
                self.l_bins_shear = self.l_shear
                self.lbin_edges_shear = None
            if self.position:
                self.l_bins_pos = self.l_pos
                self.lbin_edges_pos = None
            if self.ggl:
                self.l_bins_ggl = self.l_ggl
                self.lbin_edges_ggl = None
            if self.cmb_kappa:
                self.l_bins_kk = self.l_kk
                self.lbin_edges_kk = None
            if self.kappa_shear:
                self.l_bins_ke = self.l_ke
                self.lbin_edges_ke = None
            if self.kappa_pos:
                self.l_bins_kn = self.l_kn
                self.lbin_edges_kn = None

    def save_cls(self, block, out_path):

        if self.shear:
            if self.dobinning:
                block.put_double_array_1d(
                    'galaxy_shape_cl', 'l_bin_edges', self.lbin_edges_shear)
            block.put_double_array_1d(
                'galaxy_shape_cl', 'ell', self.l_bins_shear)
            block.put_double_array_1d(
                'galaxy_shape_cl', 'z_bin_edges', self.zbin_edges_shear)
            block.put_int('galaxy_shape_cl', 'nl', self.Nlbin_shear)
            block.put_int('galaxy_shape_cl', 'nz', self.Nzbin_shear)
        if self.position:
            if self.dobinning:
                block.put_double_array_1d(
                    'galaxy_position_cl', 'l_bin_edges', self.lbin_edges_pos)
            block.put_double_array_1d(
                'galaxy_position_cl', 'ell', self.l_bins_pos)
            block.put_double_array_1d(
                'galaxy_position_cl', 'z_bin_edges', self.zbin_edges_pos)
            block.put_int('galaxy_position_cl', 'nl', self.Nlbin_pos)
            block.put_int('galaxy_position_cl', 'nz', self.Nzbin_pos)
        if self.ggl:
            if self.dobinning:
                block.put_double_array_1d(
                    'galaxy_position_shape_cross_cl', 'l_bin_edges', self.lbin_edges_ggl)
            block.put_double_array_1d(
                'galaxy_position_shape_cross_cl', 'ell', self.l_bins_ggl)
            block.put_double_array_1d(
                'galaxy_position_shape_cross_cl', 'z_bin_edges_shear', self.zbin_edges_shear)
            block.put_double_array_1d(
                'galaxy_position_shape_cross_cl', 'z_bin_edges_position', self.zbin_edges_pos)
            block.put_int('galaxy_position_shape_cross_cl',
                          'nl', self.Nlbin_ggl)
            block.put_int('galaxy_position_shape_cross_cl',
                          'nz_shear', self.Nzbin_shear)
            block.put_int('galaxy_position_shape_cross_cl',
                          'nz_position', self.Nzbin_pos)

        out_C_ee = {}
        out_C_nn = {}
        out_C_ne = {}
        if self.shear:
            for i in range(1, self.Nzbin_shear + 1):
                for j in range(1, self.Nzbin_shear + 1):
                    bin = "bin_%d_%d" % (i, j)
                    block.put_double_array_1d(
                        'galaxy_shape_cl', bin, self.C_ee_binned[i - 1][j - 1])
                    if out_path is not None:
                        out_C_ee[bin] = self.C_ee_binned[i - 1][j - 1]
        if self.position:
            for i in range(1, self.Nzbin_pos + 1):
                for j in range(1, self.Nzbin_pos + 1):
                    bin = "bin_%d_%d" % (i, j)
                    block.put_double_array_1d(
                        'galaxy_position_cl', bin, self.C_nn_binned[i - 1][j - 1])
                    if out_path is not None:
                        out_C_nn[bin] = self.C_nn_binned[i - 1][j - 1]
        if self.ggl:
            for i in range(1, self.Nzbin_pos + 1):
                for j in range(1, self.Nzbin_shear + 1):
                    bin = "bin_%d_%d" % (i, j)
                    block.put_double_array_1d(
                        'galaxy_position_shape_cross_cl', bin, self.C_ne_binned[i - 1][j - 1])
                    if out_path is not None:
                        out_C_ne[bin] = self.C_ne_binned[i - 1][j - 1]

        try:
            omega_de = block['cosmological_parameters', 'omega_de']
        except:
            omega_k = block['cosmological_parameters', 'omega_k']
            omega_de = 1.0 - \
                block['cosmological_parameters', 'omega_m'] - omega_k
        if out_path is not "":
            cospar = {'omega_m': block['cosmological_parameters', 'omega_m'],
                      'omega_de': omega_de,
                      'omega_b': block['cosmological_parameters', 'omega_b'],
                      'h': block['cosmological_parameters', 'h0'],
                      'sigma_8': block['cosmological_parameters', 'sigma_8'],
                      'n_s': block['cosmological_parameters', 'n_s'],
                      'w0': block['cosmological_parameters', 'w'],
                      'wa': block['cosmological_parameters', 'wa']
                      }
            datavector = {'cosmology': cospar}

            if self.shear:
                datavector['C_ee'] = out_C_ee
                datavector['l_bin_edges_shear'] = self.lbin_edges_shear
                datavector['z_bin_edges_shear'] = self.zbin_edges_shear
                datavector['ell_shear'] = self.l_bins_shear
            if self.position:
                datavector['C_nn'] = out_C_nn
                datavector['l_bin_edges_pos'] = self.lbin_edges_pos
                datavector['z_bin_edges_pos'] = self.zbin_edges_pos
                datavector['ell_pos'] = self.l_bins_pos
            if self.ggl:
                datavector['C_ne'] = out_C_ne
                datavector['l_bin_edges_ggl'] = self.lbin_edges_ggl
                datavector['ell_ggl'] = self.l_bins_ggl
            import pickle as pickle

            f = open(out_path, 'wb')
            print('Saved datavector to %s' % out_path)
            pickle.dump(datavector, f)
            f.close()


def get_binned_cl(Cl, l, lbin_edges, dobinning, window):
    if dobinning:
        Cl_binned = np.zeros(len(lbin_edges) - 1)
        for i in range(len(lbin_edges) - 1):
            lmin = lbin_edges[i]
            lmax = lbin_edges[i + 1]
            sel = (l > lmin) & (l < lmax)
            if window == 'tophat':
                Cl_binned[i] = mstats.gmean(Cl[sel])
                if Cl_binned[i] == 0.0:
                    print("WARNING: Binned power spectrum is exactly 0.0.")
                if not np.isfinite(Cl_binned[i]):
                    print("WARNING: Binned power spectrum contains infinities.")
                    import pdb
                    pdb.set_trace()
            elif window == 'tophat-arithmetic':
                Cl_binned[i] = np.mean(Cl[sel])
            elif window == 'delta':
                i0 = len(Cl[sel]) / 2
                Cl_binned[i] = Cl[sel][i0]
        return Cl_binned
    else:
        return Cl


def get_binned_number_densities(nzbin, ngal):
    """
    Calculate the average number density of galaxies in each redshift bin,
    assuming an equal number in each. 
    """
    n_binned = []
    # Convert number density from arcmin^-2 to sr^-1
    ngal = (60 * 60 * 180 * 180 / (np.pi * np.pi)) * ngal
    for i in range(nzbin):
        n_binned += [ngal / nzbin]
    n_binned = np.array(n_binned)

    return n_binned
