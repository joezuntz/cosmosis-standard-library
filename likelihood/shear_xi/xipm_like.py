from __future__ import print_function
from builtins import range
from builtins import object
import numpy as np
import os
import scipy.interpolate as interpolate

dirname = os.path.split(__file__)[0]
DEFAULT_COVMAT = os.path.join(dirname, 'covariance_matrix.dat')
DEFAULT_DATA = os.path.join(dirname, 'cfhtlens_xipm_6bin.dat')


class XipmLikelihood(object):
    def __init__(self, covmat_file, xipm_data, n_z_bins, cov_num_rlzn=None, theta_mins=None, theta_maxs=None, plus_only=False, pmc=False):
        self.xipm_data = xipm_data
        # self.covmat_orig=covmat
        self.n_z_bins = n_z_bins
        self.n_z_bin_pairs = self.n_z_bins * \
            (self.n_z_bins + 1) / 2  # number of z bin pairs
        self.theta_mins = theta_mins
        self.theta_maxs = theta_maxs
        self.pmc = pmc
        if (self.theta_mins is not None) and (self.theta_maxs is None):
            self.theta_maxs = [1e9] * len(self.theta_mins)
        if (self.theta_maxs is not None) and (self.theta_mins is None):
            self.theta_mins = [0.] * len(self.theta_maxs)
        self.plus_only = plus_only
        self.thetas_arcm, self.data_vector, self.cut_inds = self.get_data_vector()
        self.n_theta = len(xipm_data['theta'])
        self.n_tot = self.n_z_bin_pairs * self.n_theta * 2
        self.covmat_orig = self.load_covariance(covmat_file)
        self.covmat = self.cut_covariance()

        weight = np.linalg.inv(self.covmat)
        p = len(self.data_vector)
        print('Using xi plus/minus: %d data points' % p)
        assert weight.shape[0] == weight.shape[1] == p
        if cov_num_rlzn is not None:  # Apply Anderson-Hartlap correction to inverse covariance
            n_mu = cov_num_rlzn
            alpha = (n_mu - p - 2.0) / (n_mu - 1.0)
            weight *= alpha
            print('applying Anderson-Hartlapp correction: ', alpha)
        self.weight = weight
        # np.savetxt('inv_cov.txt',weight)

    def load_covariance(self, filename):
        pmc = self.pmc
        try:
            cov = np.load(filename)
            return cov
        except IOError:
            try:
                data = np.loadtxt(filename).T
            except IOError:
                print("couldn't load covariance as text file or .npy array")
                return 1
        assert data.shape[0] == data.shape[1] == self.n_tot
        if pmc:
            print('Covariance in pmc format.')
            cov = np.zeros_like(data)
            k1 = 0
            for i1 in range(0, self.n_z_bins):
                for j1 in range(0, self.n_z_bins):
                    if j1 < i1:
                        continue
                    k2 = 0
                    for i2 in range(0, self.n_z_bins):
                        for j2 in range(0, self.n_z_bins):
                            if j2 < i2:
                                continue
                            cov_pp = data[2 * k1 * self.n_theta:(
                                2 * k1 + 1) * self.n_theta, 2 * k2 * self.n_theta:(2 * k2 + 1) * self.n_theta]
                            # print cov_pp
                            cov_mm = data[(2 * k1 + 1) * self.n_theta:(2 * k1 + 2) * self.n_theta,
                                          (2 * k2 + 1) * self.n_theta:(2 * k2 + 2) * self.n_theta]
                            cov_pm = data[2 * k1 * self.n_theta:(
                                2 * k1 + 1) * self.n_theta, (2 * k2 + 1) * self.n_theta:(2 * k2 + 2) * self.n_theta]
                            cov_mp = data[(2 * k1 + 1) * self.n_theta:(2 * k1 + 2) *
                                          self.n_theta, 2 * k2 * self.n_theta:(2 * k2 + 1) * self.n_theta]
                            cov[k1 * self.n_theta:(k1 + 1) * self.n_theta, k2 *
                                self.n_theta:(k2 + 1) * self.n_theta] = cov_pp
                            cov[(self.n_z_bin_pairs + k1) * self.n_theta:(self.n_z_bin_pairs + k1 + 1) * self.n_theta,
                                (self.n_z_bin_pairs + k2) * self.n_theta:(self.n_z_bin_pairs + k2 + 1) * self.n_theta] = cov_mm
                            cov[k1 * self.n_theta:(k1 + 1) * self.n_theta, (self.n_z_bin_pairs + k2) * self.n_theta:(
                                self.n_z_bin_pairs + k2 + 1) * self.n_theta] = cov_pm
                            cov[(self.n_z_bin_pairs + k1) * self.n_theta:(self.n_z_bin_pairs + k1 + 1)
                                * self.n_theta, k2 * self.n_theta:(k2 + 1) * self.n_theta] = cov_mp
                            k2 += 1
                    k1 += 1
            return cov
        return data

    def get_data_vector(self):
        # Assumes you have all cross correlations, saved in a dictionary/structured array or something....
        # With columns 'theta', 'xip_1_1', 'xip_1_2',...,'xim_1_1','xim_1_2' etc.
        # theta_min is a minimum theta to use for each column
        theta = self.xipm_data['theta']
        ntheta = len(theta)
        thetas = {}  # Dictionary of theta values for each bin combination
        k = 0
        cut_inds_full_p = []
        cut_inds_full_m = []
        data_vector_p = []
        data_vector_m = []
        thetas_full_p = []
        thetas_full_m = []
        for i in range(1, self.n_z_bins + 1):
            for j in range(i, self.n_z_bins + 1):
                bin_comb = '%d_%d' % (j, i)
                if self.plus_only:
                    if self.theta_mins is not None or self.theta_maxs is not None:
                        theta_min_p = self.theta_mins[k]
                        theta_max_p = self.theta_maxs[k]
                        cut_inds_p = (
                            np.where((theta < theta_min_p) | (theta > theta_max_p)))[0]
                        keep_inds_p = (
                            np.where((theta >= theta_min_p) & (theta <= theta_max_p)))[0]
                        thetas[(j, i)] = [theta[keep_inds_p], []]
                        thetas_full_p += list(theta[keep_inds_p])
                        data_vector_p += list(
                            self.xipm_data['xip_' + bin_comb])
                        cut_inds_p += k * ntheta
                        keep_inds_p += k * ntheta
                        cut_inds_full_p += list(cut_inds_p)
                    else:
                        thetas[(j, i)] = [theta, []]
                        data_vector_p += list(
                            self.xipm_data['xip_' + bin_comb])

                else:
                    if (self.theta_mins is not None) or (self.theta_maxs is not None):
                        theta_min_p, theta_min_m = self.theta_mins[k], self.theta_mins[k +
                                                                                       self.n_z_bin_pairs]
                        theta_max_p, theta_max_m = self.theta_maxs[k], self.theta_maxs[k +
                                                                                       self.n_z_bin_pairs]
                        cut_inds_p, cut_inds_m = ((np.where((theta < theta_min_p) | (theta > theta_max_p)))[0],
                                                  (np.where((theta < theta_min_m) | (theta > theta_max_m)))[0])
                        keep_inds_p, keep_inds_m = ((np.where((theta >= theta_min_p) & (theta <= theta_max_p)))[0],
                                                    (np.where((theta >= theta_min_m) & (theta <= theta_max_m)))[0])
                        # print 'cut_inds_p,cut_inds_m',cut_inds_p,cut_inds_m
                        # print 'keep_inds_p, keep_inds_m',keep_inds_p, keep_inds_m
                        thetas[(j, i)] = [
                            theta[keep_inds_p], theta[keep_inds_m]]
                        thetas_full_p += list(theta[keep_inds_p])
                        thetas_full_m += list(theta[keep_inds_m])
                        data_vector_p += list(
                            self.xipm_data['xip_' + bin_comb])
                        data_vector_m += list(
                            self.xipm_data['xim_' + bin_comb])
                        cut_inds_p += k * ntheta
                        keep_inds_p += k * ntheta
                        cut_inds_m += (self.n_z_bin_pairs + k) * ntheta
                        keep_inds_m += (self.n_z_bin_pairs + k) * ntheta
                        cut_inds_full_p += list(cut_inds_p)
                        cut_inds_full_m += list(cut_inds_m)
                    else:
                        thetas[(j, i)] = [theta, theta]
                        data_vector_p += list(
                            self.xipm_data['xip_' + bin_comb])
                        data_vector_m += list(
                            self.xipm_data['xim_' + bin_comb])

                k += 1

        cut_inds_full_p, cut_inds_full_m = np.array(
            cut_inds_full_p), np.array(cut_inds_full_m)
        data_vector_p, data_vector_m = np.array(
            data_vector_p), np.array(data_vector_m)
        cut_inds_full, data_vector = np.concatenate(
            (cut_inds_full_p, cut_inds_full_m)), np.concatenate((data_vector_p, data_vector_m))
        # print 'cut_inds_full',cut_inds_full
        cut_data_vector = np.delete(data_vector, cut_inds_full)
        return thetas, cut_data_vector, cut_inds_full

    def cut_covariance(self):
        if len(self.cut_inds) > 0:
            cov = np.delete(self.covmat_orig, self.cut_inds, axis=0)
            cov = np.delete(cov, self.cut_inds, axis=1)
            return cov
        else:
            return self.covmat_orig

    def interpolate_to_bin(self, theta, xi):
        return np.interp(self.theta, theta, xi)

    def __call__(self, xi_theory_dict):

        # xi_data is a dictionary containing all the bin pair data vectors
        xi_plus_vector = []
        xi_minus_vector = []
        # loop through the bins loading the theory data
        for i in range(1, self.n_z_bins + 1):
            for j in range(i, self.n_z_bins + 1):
                # print 'loading theory for bin (%d,%d)' % (j,i)
                # try both orderings for flexibility
                try:
                    theta_theory, xi_theory = xi_theory_dict[(i, j)]
                except KeyError:
                    theta_theory, xi_theory = xi_theory_dict[(j, i)]
                theta_theory = theta_theory[:len(theta_theory) / 2]
                xi_plus_theory = xi_theory[:len(theta_theory)]
                xi_minus_theory = xi_theory[len(theta_theory):]
                # interpolate to data values - this is correct if self.theta values are mean for each bin
                # Assume theory angles in radians
                theta_theory_arcm = np.degrees(theta_theory) * 60.
                theta_data_p, theta_data_m = self.thetas_arcm[(
                    j, i)][0], self.thetas_arcm[(j, i)][1]
                xi_plus_binned = np.interp(
                    theta_data_p, theta_theory_arcm, xi_plus_theory)
                xi_minus_binned = np.interp(
                    theta_data_m, theta_theory_arcm, xi_minus_theory)
                # print  xi_plus_binned,xi_minus_binned
                # print 'interpolating theory xip at angles: ',theta_data_p
                # print 'interpolating theory xim at angles: ',theta_data_m
                # build up longer vector of data
                xi_plus_vector.append(xi_plus_binned)
                xi_minus_vector.append(xi_minus_binned)
        # flatten to single vector
        xi_plus_vector = np.concatenate(xi_plus_vector)
        xi_minus_vector = np.concatenate(xi_minus_vector)
        # concatenate xi_plus_vector and xi_minus_vector to (hopefully) match form of data vector
        if self.plus_only:
            xi_vector = xi_plus_vector
        else:
            xi_vector = np.concatenate((xi_plus_vector, xi_minus_vector))
        assert xi_vector.shape == self.data_vector.shape
        # get chi2 and return log like

        delta = xi_vector - self.data_vector
        # print delta
        chi2 = np.dot(delta, np.dot(self.weight, delta))
        r = np.random.randn(xi_vector.size)
        sim = xi_vector + np.dot(self.covmat, r)
        return -chi2 / 2.0, xi_vector, self.data_vector, self.weight, self.covmat, sim


# ordering is (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) (2,2) (2,3) ... (5,5), (5,6), (6,6)
