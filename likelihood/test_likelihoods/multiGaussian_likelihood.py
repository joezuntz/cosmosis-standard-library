"""
MultiGaussian_likelihood.py: a test problem for our samplers.
We sample a 10 dimensional gaussian.
"""
import os
import sys
import numpy as np
import os.path
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


class Multigaussian(object):

    def __init__(self, ndim, means=None, covariances=None):
        self.ndim = ndim
        if means:
            self.means = np.loadtxt(means)
            # print self.means
            # print "loading saved means"
        else:
            self.means = np.random.rand(ndim)
            np.savetxt('./demos/means10D.out', self.means)

        if covariances:
            cov = np.loadtxt(covariances)
            # print "loading saved covariance matrix"

        else:
            cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
            cov = np.triu(cov)
            cov += cov.T - np.diag(cov.diagonal())
            cov = np.dot(cov, cov)
            np.savetxt('./demos/covariance10D.out', cov)

        self.icov = np.linalg.inv(cov)
        self.det = np.linalg.det(cov)
        self.norm = 0.5 * np.log((2.0 * np.pi) ** self.ndim * self.det)

    def lnprob(self, x):

        diff = x - self.means

        return -np.dot(diff, np.dot(self.icov, diff)) / 2.0 - self.norm

    def __call__(self, x):
        return lnprob(self, x)


def setup(options):

    try:
        means = options[option_section, "means_filename"]
        print "Existing means (%s) will be used" % means
    except Exception, e:
        means = None
        print "new means will be generated"

    try:
        covariances = options[option_section, "cov_filename"]
        print "Existing covariances (%s) will be used" % covariances
    except Exception, e:
        covariances = None
        print "new covariance matrix will be generated"

    return means, covariances


def execute(block, config):
    x = []
    x.append(block[cosmo, "MG1"])
    x.append(block[cosmo, "MG2"])
    x.append(block[cosmo, "MG3"])
    x.append(block[cosmo, "MG4"])
    x.append(block[cosmo, "MG5"])
    x.append(block[cosmo, "MG6"])
    x.append(block[cosmo, "MG7"])
    x.append(block[cosmo, "MG8"])
    x.append(block[cosmo, "MG9"])
    x.append(block[cosmo, "MG10"])
    mg = Multigaussian(len(x), config[0], config[1])
    # compute likelihood
    like = mg.lnprob(x)
    block[likes, "MULTIGAUSSIAN_LIKE"] = like

    # except KeyboardInterrupt:
    # If the user presses Ctrl-C we respect that
    #     raise KeyboardInterrupt
    # except Exception as E:
    # But any other kind of error is a problem
    #     print "There was an error calculating the MultiGaussian likelihood: "
    #     return 1

    return 0


def cleanup(config):
# nothing to do here!  We just include this
# for completeness
    return 0
