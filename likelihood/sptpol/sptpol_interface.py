from __future__ import print_function
from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
import os
import numpy as np
import sptpol_like

ROOT_DIR = os.path.split(os.path.abspath(__file__))[0]
DATA_DIR = os.path.join(ROOT_DIR, "data")
DATASET = "sptpol_500d_TEEE"
DATASET_DIR = os.path.join(DATA_DIR, DATASET)
DEFAULT_BANDPOWERS = os.path.join(DATASET_DIR, "sptpol_500d_TTTEEE.bp_file")
DEFAULT_COVMAT = os.path.join(DATASET_DIR, "sptpol_500d_TEEE.cov_file")
DEFAULT_BEAMERRS = os.path.join(DATASET_DIR, "sptpol_500d_TEEE_beamerrs.txt")
DEFAULT_WINDOWS = os.path.join(DATASET_DIR, "bpwfs_ascii")


class SPTPolLikelihood(GaussianLikelihood):
    like_name = "sptpol"

    def __init__(self, options):
        self.use_te = options.get_bool("use_te", default=True)
        self.use_ee = options.get_bool("use_ee", default=True)

        data_filename = options.get_string("bandpowers", DEFAULT_BANDPOWERS)
        cov_filename = options.get_string("covmat", DEFAULT_COVMAT)
        beamerr_file = options.get_string("beam_errors", DEFAULT_BEAMERRS)
        windows_dir = options.get_string("windows", DEFAULT_WINDOWS)

        if (not os.path.exists(DATASET_DIR)) and (cov_filename == DEFAULT_COVMAT or windows_dir == DEFAULT_WINDOWS):
            raise ValueError(
                "Please see the file {}/readme.txt for details on downloading the SPTpol data ({} does not exist)".format(ROOT_DIR))

        self.spt_data = sptpol_like.SPTPolData(
            data_filename, cov_filename, self.use_te, self.use_ee)
        self.spt_model = sptpol_like.SPTPolTheoryModel(
            beamerr_file, windows_dir, self.use_te, self.use_ee)

        print("Will look for these nuisance parameters in a section [sptpol]:")
        print("    " + (', '.join(self.spt_model.nuisance_parameter_names())))

        super(SPTPolLikelihood, self).__init__(options)

    def build_data(self):
        return None, self.spt_data.vector

    def build_covariance(self):
        return self.spt_data.covmat

    def extract_theory_points(self, block):
        # just get all supplied nuisance parameters
        nuisance = {}
        for name in self.spt_model.nuisance_parameter_names():
            nuisance[name] = block["sptpol", name]

        ell = block[names.cmb_cl, "ell"]
        te = block[names.cmb_cl, "te"]
        ee = block[names.cmb_cl, "ee"]

        # The likelihood code I wrote assumes that the spectra are all
        # suck that C_ell = C[ell].  This is only true if ell goes down to zero
        # also need to check that ell == arange(n_ell)
        if ell.min() != 0:
            start_ell = np.arange(ell.min(), dtype=np.int32)
            start_cl = np.zeros(len(start_ell))
            ell = np.concatenate([start_ell, ell])
            te = np.concatenate([start_cl, te])
            ee = np.concatenate([start_cl, ee])

        if ell.max() < self.spt_model.ell.max() + 1:
            raise ValueError("""
You need to calculate the CMB spectra to ell>{} to use the SPTpol data.
This setting can be changed in your camb/class section.""".format(self.spt_model.ell.max() + 1))
        if not np.all(ell == np.arange(len(ell))):
            raise ValueError("""
For some reason your ell values do not include all the integers - there are gaps.
The default modules should not do this so probably you messed with something. Raise an issue if not.""")

        # and the CMB
        cmb = {
            "ell": block[names.cmb_cl, "ell"],
            "te": block[names.cmb_cl, "te"],
            "ee": block[names.cmb_cl, "ee"],
        }

        theory_vector, components_te, components_ee = self.spt_model(
            cmb, nuisance)

        # In case interesting, save the different components to the block
        if self.use_te:
            for key, val in list(components_te.items()):
                if key != "ell":
                    key = key + "_te"
                block['spt_model', key] = val

        return theory_vector


setup, execute, cleanup = SPTPolLikelihood.build_module()
