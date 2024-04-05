import os
import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.gaussian_likelihood import GaussianLikelihood

ROOT_DIR = os.path.split(os.path.abspath(__file__))[0]

class DESIY1BAOLikelihood(GaussianLikelihood):

    like_name = "desi_y1_bao"

    def __init__(self, options):

        iso_file = os.path.join(ROOT_DIR, "desi_y1_bao_iso.txt")
        ani_file = os.path.join(ROOT_DIR, "desi_y1_bao_ani.txt")

        iso_cols = ['z_eff', 'DVrd', 'DVrd_error']
        ani_cols = ['z_eff', 'DMrd', 'DMrd_error', 'DHrd', 'DHrd_error', 'corr']

        self.iso_data = {k: v for k,v in zip(iso_cols, np.loadtxt(iso_file, unpack=True))}
        self.ani_data = {k: v for k,v in zip(ani_cols, np.loadtxt(ani_file, unpack=True))}

        self.niso = len(self.iso_data['z_eff'])
        self.nani = len(self.ani_data['z_eff'])

        self.feedback = options.get_bool("feedback", default=False)

        super().__init__(options)

    def build_data(self):
        zeff = np.concatenate([self.iso_data['z_eff'], self.ani_data['z_eff'], self.ani_data['z_eff']])
        dist = np.concatenate([self.iso_data['DVrd'],  self.ani_data['DMrd'],  self.ani_data['DHrd']])

        if self.feedback:
            print('Data distances')
            print('zeff   DV/rd   DM/rd   DH/rd')
            for z, v in zip(self.iso_data['z_eff'], self.iso_data['DVrd']):
                print(f'{z:.2f}   {v:.2f}')
            for z, m, h in zip(self.ani_data['z_eff'], self.ani_data['DMrd'], self.ani_data['DHrd']):
                print(f'{z:.2f}           {m:.2f}   {h:.2f}')

        return zeff, dist

    def build_covariance(self):
        cov = np.diag(np.concatenate([self.iso_data['DVrd_error'], self.ani_data['DMrd_error'], self.ani_data['DHrd_error']]))**2

        zipped = zip(self.ani_data['DMrd_error'], self.ani_data['DHrd_error'], self.ani_data['corr'])

        for i, (DMrd_error, DHrd_error, corr) in enumerate(zipped, start=self.niso):
            cov[i,i+self.nani] = DMrd_error * DHrd_error * corr
            cov[i+self.nani,i] = cov[i,i+self.nani]

        return cov

    def build_inverse_covariance(self):
        return np.linalg.inv(self.cov)

    def extract_theory_points(self, block):

        z = block[names.distances, 'z']

        # Sound horizon at the drag epoch
        rd = block[names.distances, "rs_zdrag"]
        if self.feedback:
            print(f'rs_zdrag = {rd}')

        # Comoving distance
        DM_z = block[names.distances, 'd_m']  # in Mpc

        # Hubble distance
        DH_z = 1/block[names.distances, 'H'] # in Mpc

        # Angle-averaged distance
        DV_z = (z * DM_z**2 * DH_z)**(1/3) # in Mpc

        # z and distance maybe are loaded in chronological order
        # Reverse to start from low z
        if (z[1] < z[0]):
            z  = z[::-1]
            DM_z = DM_z[::-1]
            DH_z = DH_z[::-1]

        # Find theory DM and DH at effective redshift by interpolation
        z_eff = self.data_x[:self.niso+self.nani]

        DMrd = np.interp(z_eff, z, DM_z)/rd
        DHrd = np.interp(z_eff, z, DH_z)/rd
        DVrd = (z_eff * DMrd**2 * DHrd)**(1/3)

        if self.feedback:
            print('Theory distances')
            print('zeff   DV/rd   DM/rd   DH/rd')
            for i, (z, m, h, v) in enumerate(zip(z_eff, DMrd, DHrd, DVrd)):
                if i < self.niso:
                    print(f'{z:.2f}   {v:.2f}')
                else:
                    print(f'{z:.2f}           {m:.2f}   {h:.2f}')

        return np.concatenate([DVrd[:self.niso], DMrd[self.niso:], DHrd[self.niso:]])

setup, execute, cleanup = DESIY1BAOLikelihood.build_module()

