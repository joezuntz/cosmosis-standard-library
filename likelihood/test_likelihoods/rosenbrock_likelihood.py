"""
rosenbrock_likelihood.py: a test problem for our samplers.
We sample a Rosenbrock function.
"""

from builtins import object
import os
import sys
import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


class Rosenbrock(object):

    def __init__(self, a1=100., a2=20.):
        self.a1 = a1
        self.a2 = a2

    def loglikelihood(self, p):
        return -(self.a1 * (p[1] - p[0] ** 2) ** 2 + (1 - p[0]) ** 2) / self.a2

    def __call__(self, p):
        return self.loglikelihood(p)


def setup(options):
    a1 = options.get_double(option_section, "a1", default=100.)
    a2 = options.get_double(option_section, "a2", default=20.)
    return Rosenbrock(a1, a2)


def execute(block, rosenbrock):
    p = [block[cosmo, "P1"], block[cosmo, "P2"]]
    block[likes, 'ROSENBROCK_LIKE'] = rosenbrock(p)
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    pass
