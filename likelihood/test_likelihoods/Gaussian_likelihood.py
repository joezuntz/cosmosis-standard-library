from __future__ import print_function
from builtins import object
import os
import pydesglue
import numpy as np
try:
    import pylab
except ImportError:
    pylab = None


class Multigaussian(object):
    def __init__(self, ndim):
        self.ndim = ndim

        self.means = 0.0  # *np.random.rand(ndim)
        self.cov = 0.5  # -np.random.rand(ndim**2).reshape((ndim, ndim)
        self.norm = 0.5 * np.log((2.0 * np.pi)**self.ndim * self.cov)

    def lnprob(self, x):

        diff = x[0] - self.means
        print("x", x[0], diff**2, type(self.cov))

        return -(diff**2) / (self.cov * 2.0) - self.norm

    def __call__(self, x):
        return lnprob(self, x)


def execute(handle):
    try:
        package = pydesglue.DesDataPackage.from_fits_handle(handle)
        section = pydesglue.section_names.cosmological_parameters
        x = []
        x.append(package.get_param(section, "MG1"))
        mg = Multigaussian(len(x))
        print("here")
        # compute likelihood
        like = mg.lnprob(x)
        section = pydesglue.section_names.likelihoods
        package.set_param(section, "GAUSSIAN_LIKE", like)
        package.write_to_fits_handle(handle)

    except KeyboardInterrupt:
        # If the user presses Ctrl-C we respect that
        raise KeyboardInterrupt
    except Exception as E:
        # But any other kind of error is a problem
        print("There was an error calculating the MultiGaussian likelihood: ")
        return 1
        return 0
