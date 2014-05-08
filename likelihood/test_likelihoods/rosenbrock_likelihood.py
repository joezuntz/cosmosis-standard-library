"""
rosenbrock_likelihood.py: a test problem for our samplers.
We sample a Rosenbrock function.
"""

import os
import sys
import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


class Rosenbrock(object):

    def __init__(self):
        self.a1 = 100.0
        self.a2 = 20.0

    def logl(self, p):
        return -(self.a1 * (p[1] - p[0] ** 2) ** 2 + (1 - p[0]) ** 2) / self.a2

    def __call__(self, p):
        return np.log(self, p)


def setup(options):
    return 1


def execute(block, config):
    p = []
    p.append(block[cosmo, "P1"])
    p.append(block[cosmo, "P2"])

    rosenbrock = Rosenbrock()
    # compute likelihood
    like = rosenbrock.logl(p)
    block[likes, 'ROSENBROCK_LIKE'] = like

    # except KeyboardInterrupt :
    # If the user presses Ctrl-C we respect that
    #     raise

    # except Exception as E:
    # But any other kind of error is a problem
    #     print "There was an error calculating the Rosenbrock likelihood: "
    #     return 1

    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
