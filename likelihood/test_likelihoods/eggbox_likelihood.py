"""
eggbox_likelihood.py: a test problem for our samplers.
We sample an eggbox function a very pathological case.
"""

from builtins import object
import os
import sys
import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods


class Eggbox(object):

    """Adapted from https://github.com/dfm/emcee/blob/master/examples/eggbox.py """

    def __init__(self, tmax=10. * np.pi):
        self.tmax = tmax

    def loglhood(self, t):
        return (2.0 + np.cos(t[0] / 2.0) * np.cos(t[1] / 2.0)) ** 5.0

    def __call__(self, t):
        return self.loglhood(t)


def setup(options):
    tmax = options.get_double(option_section, "tmax", default=10. * np.pi)
    return Eggbox(tmax)


def execute(block, eggbox):
    t = [block[cosmo, 'T1'],
         block[cosmo, 'T2']]
    block[likes, 'EGGBOX_LIKE'] = eggbox(t)
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
