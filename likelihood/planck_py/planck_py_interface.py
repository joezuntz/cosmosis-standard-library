from cosmosis.datablock import names
from cosmosis.gaussian_likelihood import GaussianLikelihood
import planck_lite_py
import numpy as np
import os

dirname = os.path.split(__file__)[0]

class PlanckPythonLikelihood(GaussianLikelihood):
    x_section = 'cmb_cl'
    x_name = 'ell'
    y_section = 'cmb_cl'
    like_name = 'planck'

    def build_data(self):
        year = self.options['year']
        spectra = self.options['spectra']
        use_low_ell_bins = self.options['use_low_ell_bins']

        # Use the local directory by default
        data_directory = os.path.join(dirname, 'data')
        data_directory = self.options.get_string('data_directory', data_directory)

        self.calculator = planck_lite_py.PlanckLitePy(year=year, spectra=spectra, 
            use_low_ell_bins=use_low_ell_bins, data_directory=data_directory)
        x = None
        y = self.calculator.mu
        return x, y

    def build_covariance(self):
        return self.calculator.cov

    def build_inverse_covariance(self):
        return self.calculator.fisher

    def extract_theory_points(self, block):
        ell = block[self.x_section, self.x_name]
        ellmin = ell[0]

        Dltt = block[self.y_section, 'TT']
        Dlte = block[self.y_section, 'TE']
        Dlee = block[self.y_section, 'EE']

        y = self.calculator.make_mean_vector(Dltt, Dlte, Dlee, ellmin=ellmin)
        y = self.calculator._cut_vector(y)

        return y

setup, execute, cleanup = PlanckPythonLikelihood.build_module()
