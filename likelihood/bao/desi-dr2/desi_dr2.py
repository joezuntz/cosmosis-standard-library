import os
from cosmosis.gaussian_likelihood import GaussianLikelihood
import numpy as np
import scipy.interpolate
import scipy.linalg
from astropy.table import Table

# There are now only two kinds of measurement here.
# All the measurements except one now include measurements
# of DV/rd, DM/DH, DM/rd, and DH/rd. The one exception is
# the BGS measurements which is just DV/rd

# We have marked the other measurements as all NaNs
# in the data to indicate this

THIS_DIR = os.path.split(__file__)[0]
DR2_DEFAULT_FILENAME = os.path.join(THIS_DIR, "dr2_data.txt")
KIND1_DV_RD = "DV/rd"
KIND2_DM_DH = "DM/DH"
KIND3_DM_RD = "DM/rd"
KIND4_DH_RD = "DH/rd"


def extract_desi_table_data_vector_covariance(data):
    z = []
    tracers = []
    kinds = []
    data_vector = []
    cov_blocks = []

    for row in data:
        tracer = row["Tracer"]
        zeff = row["zeff"]

        if np.isnan(row["DM/DH"]):
            # The BGS sample only measures DV/r_d, and so we
            # assume that as in the default data file this is labelled
            # with NaNs in the other columns
            z.append(zeff)
            tracers.append(tracer)
            kinds.append(KIND1_DV_RD)
            data_vector.append(row["DV/rd"])
            c = np.zeros((1, 1))
            c[0, 0] = row["DV/rd_err"]**2
            cov_blocks.append(c)
        else:
            z.extend([zeff, zeff, zeff, zeff])
            tracers.extend([tracer, tracer, tracer, tracer])
            kinds.extend([KIND1_DV_RD, KIND2_DM_DH, KIND3_DM_RD, KIND4_DH_RD])
            data_vector.extend([row[KIND1_DV_RD], row[KIND2_DM_DH], row[KIND3_DM_RD], row[KIND4_DH_RD]])
            c12 = row["rV,M/H"]
            c34 = row["rM,H"]
            sigma1 = row["DV/rd_err"]
            sigma2 = row["DM/DH_err"]
            sigma3 = row["DM/rd_err"]
            sigma4 = row["DH/rd_err"]
            c = np.zeros((4, 4))
            c[0, 0] = sigma1**2
            c[1, 1] = sigma2**2
            c[2, 2] = sigma3**2
            c[3, 3] = sigma4**2
            c[0, 1] = c[1, 0] = c12*sigma1*sigma2
            c[2, 3] = c[3, 2] = c34*sigma3*sigma4
            cov_blocks.append(c)
    
    covariance = scipy.linalg.block_diag(*cov_blocks)
    return np.array(z), np.array(tracers), np.array(kinds), np.array(data_vector), covariance


class DESILikelihood(GaussianLikelihood):
    """
    The 2025 DESI DR2 likelihoods from https://arxiv.org/pdf/2503.14738

    We allow the user to specify which data sets to use, and combine
    them all into one. The data sets are:
    - BGS
    - LRG1
    - LRG2
    - LRG3+ELG1
    - ELG2
    - QSO
    - Lya QSO

    Two data sets are available but not included when using the "all"
    setting because they are correlated with and superseded by the
    combined version:

    - LRG3
    - ELG1

    """
    # users can override this if they want to use a different name
    # which can be useful if you want to keep the different likelihoods
    # separately.
    like_name = "desi_bao"
    x_section = 'distances'
    x_name = 'z'
    y_section = 'distances'

    def __init__(self, options):
        data_file = options.get_string("data_file", DR2_DEFAULT_FILENAME)
        self.table = Table.read(data_file, format="ascii.commented_header")

        data_sets = options.get_string("desi_data_sets", default="all")
        available_tracers = list(self.table["Tracer"])

        if data_sets == "all":
            data_sets = available_tracers[:]
            data_sets.remove("LRG3")
            data_sets.remove("ELG1")
        else:
            data_sets = data_sets.replace(",", " ").split()

        for data_set in data_sets:
            data_set = data_set.strip()
            if data_set not in available_tracers:
                raise ValueError(f"Unknown DESI data set {data_set}. Valid options are: {available_tracers} (space or comma-separated to use more than one)")

        rows = [available_tracers.index(data_set) for data_set in data_sets]
        self.table = self.table[rows]
        print("Using these DESI data sets:\n", self.table["Tracer"])

        super().__init__(options)
    

    def build_data(self):
        z, tracers, kinds, mu, C = extract_desi_table_data_vector_covariance(self.table)
        self.kinds = kinds
        self.tracers = tracers
        self.table_cov = C
        

        # record the indices of the d_v and d_m_d_h measurements
        # for later
        self.dv_rd_index = np.where(kinds==KIND1_DV_RD)[0]
        self.dm_dh_index = np.where(kinds==KIND2_DM_DH)[0]
        self.dm_rd_index = np.where(kinds==KIND3_DM_RD)[0]
        self.dh_rd_index = np.where(kinds==KIND4_DH_RD)[0]

        self.need_dv = len(self.dv_rd_index) > 0
        self.need_dm = len(self.dm_dh_index) > 0 or len(self.dm_rd_index) > 0
        self.need_dh = len(self.dh_rd_index) > 0 or len(self.dm_dh_index) > 0

        return z, mu


    def build_covariance(self):
        # This has already been calculated in build_data
        return self.table_cov

    def extract_theory_points(self, block):
        z_theory = block[self.x_section, self.x_name]
        y = np.zeros(self.data_x.size)
        r_s = block[self.y_section, "rs_zdrag"]
        if self.need_dv:
            d_v = block[self.y_section, "d_v"]
        if self.need_dm:
            d_m = block[self.y_section, "d_m"]
        if self.need_dh:
            d_h = 1.0 / block[self.y_section, "h"]

        block["distances", "h0rd"] = block["cosmological_parameters", "h0"] * r_s

        if len(self.dv_rd_index):
            z_data = self.data_x[self.dv_rd_index]
            f = scipy.interpolate.interp1d(z_theory, d_v / r_s, kind=self.kind)
            y[self.dv_rd_index] = f(z_data)

        if len(self.dm_dh_index):
            z_data = self.data_x[self.dm_dh_index]
            f = scipy.interpolate.interp1d(z_theory, d_m / d_h, kind=self.kind)
            y[self.dm_dh_index] = f(z_data)

        if len(self.dm_rd_index):
            z_data = self.data_x[self.dm_rd_index]
            f = scipy.interpolate.interp1d(z_theory, d_m / r_s, kind=self.kind)
            y[self.dm_rd_index] = f(z_data)

        if len(self.dh_rd_index):
            z_data = self.data_x[self.dh_rd_index]
            f = scipy.interpolate.interp1d(z_theory, d_h / r_s, kind=self.kind)
            y[self.dh_rd_index] = f(z_data)

        return y

setup, execute, cleanup = DESILikelihood.build_module()
