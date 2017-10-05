from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
from scipy.interpolate import interp1d
import os
import numpy as np

dirname = os.path.split(__file__)[0]
default_data_file = os.path.join(dirname, "data.txt")
default_covmat_file = os.path.join(dirname, "covmat.npy")


class PlanckSimpleLikelihood(GaussianLikelihood):
    like_name = "planck2015_simple"
    x_section = "cmb_cl"
    y_section = "cmb_cl"
    x_name = "ell"

    def build_data(self):
        data_file = self.options.get_string("data_file", default_data_file)
        ell, c_ell = np.loadtxt(data_file).T
        te_start, ee_start = np.where(np.diff(ell) < 0)[0] + 1
        ell = np.split(ell, [te_start, ee_start])
        return ell, c_ell

    def build_covariance(self):
        # covmat should be stored in npy format for size
        covmat_file = self.options.get_string(
            "covmat_file", default_covmat_file)
        covmat = np.load(covmat_file)
        return covmat

    def extract_theory_points(self, block):
        # The data is supplied as c_ell, not as l**2 c_ell,
        # so convert our theory predictions
        ell_theory = block[self.x_section, self.x_name]
        f = ell_theory * (ell_theory + 1) / (2 * np.pi)
        tt_theory = block[self.y_section, "tt"] / f
        te_theory = block[self.y_section, "te"] / f
        ee_theory = block[self.y_section, "ee"] / f

        # interpolate into the theory,
        tt_predicted = interp1d(ell_theory, tt_theory)(self.data_x[0])
        te_predicted = interp1d(ell_theory, te_theory)(self.data_x[1])
        ee_predicted = interp1d(ell_theory, ee_theory)(self.data_x[2])

        cl_predicted = np.concatenate(
            [tt_predicted, te_predicted, ee_predicted])

        return cl_predicted


setup, execute, cleanup = PlanckSimpleLikelihood.build_module()
