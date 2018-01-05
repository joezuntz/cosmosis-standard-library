"""
MultiGaussian_likelihood.py: a test problem for our samplers.
We sample a 10 dimensional gaussian.
"""
from __future__ import print_function
from builtins import range
from builtins import object
import os
import sys
import numpy as np
import os.path
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.datablock import BlockError


params = section_names.test_parameters
likes = section_names.likelihoods


class MultiGaussian(object):
    NDIM = 10

    def __init__(self, means=None, covariances=None):
        if means:
            self.means = np.loadtxt(means)
            if len(self.means) != self.NDIM:
                raise ValueError("The number of gaussian means in %s does "
                                 "not match those in multiGaussian_likelihood.py")
        else:
            self.means = np.random.rand(self.NDIM)
            np.savetxt("demos/means%dD.out" % (self.NDIM,), self.means)

        if covariances.lower() == "identity":
            cov = np.identity(self.NDIM) * 0.1
        elif covariances is None:
            cov = 0.5 - np.random.rand(self.NDIM **
                                       2).reshape((self.NDIM, self.NDIM))
            cov = np.triu(cov)
            cov += cov.T - np.diag(cov.diagonal())
            cov = np.dot(cov, cov)
            np.savetxt('./demos/covariance10D.out', cov)
        else:
            cov = np.loadtxt(covariances)
            if cov.shape != (self.NDIM, self.NDIM):
                raise ValueError("The shape of the gaussian covariance file %s does "
                                 "not match what was expected (%d, %d) vs (%d, %d)" %
                                 (cov.shape[0], cov.shape[1], self.NDIM, self.NDIM))

        self.icov = np.linalg.inv(cov)
        self.det = np.linalg.det(cov)
        self.norm = 0.5 * np.log((2.0 * np.pi)**self.NDIM * self.det)

    def lnprob(self, x):
        diff = x - self.means
        return -np.dot(diff, np.dot(self.icov, diff)) / 2.0 - self.norm

    def __call__(self, x):
        return lnprob(self, x)


def setup(options):
    try:
        means = options[option_section, "means_filename"]
        print("Existing means (%s) will be used" % means)
    except BlockError:
        means = None
        print("new means will be generated")

    try:
        covariances = options[option_section, "cov_filename"]
        print("Existing covariances (%s) will be used" % covariances)
    except BlockError:
        covariances = None
        print("new covariance matrix will be generated")

    return MultiGaussian(means, covariances)


def execute(block, mg):
    x = [block[params, "MG%d" % (i,)] for i in range(1, mg.NDIM + 1)]
    block[likes, "MULTIGAUSSIAN_LIKE"] = mg.lnprob(x)
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
