"""
MultiGaussian_likelihood.py: a test problem for our samplers.
We sample a 10 dimensional gaussian.
"""
import os
import numpy as np
import os.path
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


class Multigaussian(object):

    def __init__(self, ndim):
        self.ndim = ndim
        if os.path.isfile('./means10D.out'):
            self.means = np.loadtxt('./means10D.out')
            # print "loading saved means"
        # else:
        #     self.means = np.random.rand(ndim)
        #     np.savetxt('./means10D.out', self.means)

        if os.path.isfile('./covariance10D.out'):
            cov = np.loadtxt('./covariance10D.out')
            # print "loading saved covariance matrix"

        # else:
        #     cov = 0.5-np.random.rand(ndim**2).reshape((ndim, ndim))
        #     cov = np.triu(cov)
        #     cov += cov.T - np.diag(cov.diagonal())
        #     cov = np.dot(cov,cov)
        #     np.savetxt('covariance10D.out', cov)

        self.icov = np.linalg.inv(cov)
        self.det = np.linalg.det(cov)
        self.norm = 0.5 * np.log((2.0 * np.pi) ** self.ndim * self.det)

    def lnprob(self, x):

        diff = x - self.means

        return -np.dot(diff, np.dot(self.icov, diff)) / 2.0 - self.norm

    def __call__(self, x):
        return lnprob(self, x)


def setup(options):
    section = option_section
    return 1


def execute(block, config):
        x = []
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        x.append(block[cosmo, "MG1"])
        mg = Multigaussian(len(x))
        # compute likelihood
        like = mg.lnprob(x)
        block[likes, "MULTIGAUSSIAN_LIKE"] = like

    except KeyboardInterrupt:
        # If the user presses Ctrl-C we respect that
        raise KeyboardInterrupt
    except Exception as E:
        # But any other kind of error is a problem
        print "There was an error calculating the MultiGaussian likelihood: "
        return 1

        return 0

def cleanup(config):
# nothing to do here!  We just include this
# for completeness
    return 0
