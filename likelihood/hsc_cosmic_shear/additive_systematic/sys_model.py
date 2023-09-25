import os
import numpy as np
import astropy.table as astTable
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section

psfsec = "psf_systematics_parameters"  # PSF field name in data block
rescale = 1.0 / 60.0 / 180.0 * np.pi  # scaling factor from arcmin to radian

def generate_delta_xip(params1, params2, sys_cors1, sys_cors2):
    """Generates the systematics error on xip using the parameters

    Args:
        params1 (ndarray):      alpha-like systematic parameters [first]
        params2 (ndarray):      alpha-like systematic parameters [second]
        sys_cors1 (ndarray):    systematic property matrix [first component]
        sys_cors2 (ndarray):    systematic property matrix [second component]
    Returns:
        delta_xip (ndarray):    the delta xip from the PSF systematic parameters
    """

    if not isinstance(sys_cors1, np.ndarray):
        raise TypeError("The input systematic matrix1 should be ndarray.")
    if not isinstance(sys_cors2, np.ndarray):
        raise TypeError("The input systematic matrix2 should be ndarray.")
    if not isinstance(params1, np.ndarray):
        raise TypeError("The input parameter1 should be ndarray.")
    if not isinstance(params2, np.ndarray):
        raise TypeError("The input parameter2 should be ndarray.")
    # sys_cors1 and sys_cors2 are in shape of (ncor, ncor, n_theta_bin)
    # for each theta bin it is
    # [
    # [ pp1 pq1 p1 ]
    # [ pq1 qq1 q1 ]
    # [ p1  q1  1  ]
    # ]
    # params is in shape of (nzs, ncor)
    # for each source redshift
    # [
    # alpha, beta, c1
    # ]
    # [
    # alpha, beta, c2
    # ]

    ncor = sys_cors1.shape[0]
    assert sys_cors1.shape[1] == ncor, "sys_cors1.shape is wrong."
    n_theta_bin = sys_cors1.shape[-1]
    assert sys_cors2.shape == (ncor, ncor, n_theta_bin), "sys_cors2.shape is wrong"
    if len(params1.shape) == 1:
        params1 = params1[None, :]
    if len(params2.shape) == 1:
        params2 = params2[None, :]

    nzs = params1.shape[0]
    if params1.shape != (nzs, ncor):
        raise ValueError(
            "The shape of parameters 1 (%s) are not correct" % (params1.shape)
        )
    if params2.shape != (nzs, ncor):
        raise ValueError(
            "The shape of parameters 2 (%s) are not correct" % (params2.shape)
        )

    return _generate_delta_xip(
        nzs,
        n_theta_bin,
        params1,
        params2,
        sys_cors1,
        sys_cors2,
    )


def _generate_delta_xip(nzs, n_theta_bin, params1, params2, sys_cors1, sys_cors2):
    ncross = ((nzs + 1) * nzs) // 2  # number of cross-correlations
    delta_xip = np.zeros(shape=(ncross, n_theta_bin))
    zcs = 0  # id of the cross correlation (from 0 to ncross)
    for zi in range(nzs):
        for zj in range(zi, nzs):
            delta_xip[zcs] = np.dot(
                params1[zj], np.dot(params1[zi], sys_cors1)
            ) + np.dot(params2[zj], np.dot(params2[zi], sys_cors2))
            zcs += 1  # updates id
    return delta_xip

class Interp1d(object):
    def __init__(self, angle, spec, bounds_error=False):
        if np.all(spec > 0):
            self.interp_func = interp1d(
                np.log(angle),
                np.log(spec),
                bounds_error=bounds_error,
                fill_value=-np.inf,
            )
            self.interp_type = "loglog"
            self.x_func = np.log
            self.y_func = np.exp
        elif np.all(spec < 0):
            self.interp_func = interp1d(
                np.log(angle),
                np.log(-spec),
                bounds_error=bounds_error,
                fill_value=-np.inf,
            )
            self.interp_type = "minus_loglog"
            self.x_func = np.log
            self.y_func = lambda y: -np.exp(y)
        else:
            self.interp_func = interp1d(
                np.log(angle), spec, bounds_error=bounds_error, fill_value=0.0
            )
            self.interp_type = "log_ang"
            self.x_func = np.log
            self.y_func = lambda y: y

    def __call__(self, angle):
        interp_vals = self.x_func(angle)
        try:
            spec = self.y_func(self.interp_func(interp_vals))
        except ValueError:
            interp_vals[0] *= 1 + 1.0e-9
            interp_vals[-1] *= 1 - 1.0e-9
            spec = self.y_func(self.interp_func(interp_vals))
        return spec

def setup(options):
    # sclae cut on xip, in units of arcmin
    theta_min = options.get_double(option_section, "theta_min", 3.5) * rescale
    theta_max = options.get_double(option_section, "theta_max", 200.0) * rescale
    model_type = options.get_string(option_section, "model_type")

    if model_type == "systematics":
        psf_file = options.get_string(option_section, "psf_file")
        nzs = options.get_int(option_section, "nzs")
        data = np.load(psf_file)
        sys1 = data["sys1"]
        sys2 = data["sys2"]
        npar = sys1.shape[0]
        # (npar, npar, ntheta)
        r = data["r"] * rescale  # from arcmin to radian
        ntheta = len(r)
        assert sys1.shape == (npar, npar, ntheta)
        assert sys2.shape == (npar, npar, ntheta)

        tfile = options.get_string(option_section, "tfile", "")
        if len(tfile) > 0:
            tobj = np.load(tfile)
            jacobian = tobj["jacobian"]
            center = tobj["center"]
        else:
            jacobian = np.eye((npar - 1) * nzs)  # A null transform
            # -1 since the last dimension is for c
            center = 0.0
        out = {
            "model_type": model_type,
            "nzs": nzs,
            "npar": npar,
            "r": r,
            "sys1": sys1,
            "sys2": sys2,
            "theta_min": theta_min,
            "theta_max": theta_max,
            "jacobian": jacobian,
            "center": center,
        }
    else:
        raise ValueError("model_type: %s does not support" % model_type)
    return out


def execute(block, config):
    n_a = block["shear_xi_plus", "nbin_a"]
    n_b = block["shear_xi_plus", "nbin_b"]
    assert n_a == n_b
    ncross = n_a * (n_a + 1) // 2
    # cut theta
    tmin = config["theta_min"]
    tmax = config["theta_max"]
    theta = block["shear_xi_plus", "theta"]
    # theta is in units on radian
    mskT = np.greater(theta, tmin) & np.less(theta, tmax)
    theta = theta[mskT]
    block["shear_xi_plus", "theta"] = theta

    # read the theta_min, theta_max and the pca object
    if config["model_type"] == "pca":
        pca = config["model"]
        npar = config["npar"]
        # read the parameters from data block
        pars = np.array([block[psfsec, "p%d" % (i + 1)] for i in range(npar)])
        ndim = pca.r.size // ncross
        dxip = pca.itransform(pars).reshape(ncross, ndim)
        rr = pca.r[:ndim]
        dxip = Interp1d(rr, dxip, bounds_error=False)(theta)
    elif config["model_type"] == "systematics":
        nzs = config["nzs"]
        if nzs != 1:
            assert nzs == n_a
        npar = config["npar"]
        sys1 = config["sys1"]
        sys2 = config["sys2"]
        # correlation parameters (alpha, beta..)
        pars_cor = np.array(
            [
                [
                    block[psfsec, "psf_cor%d_z%d" % (i + 1, j + 1)]
                    for i in range(npar - 1)
                ]
                for j in range(nzs)
            ]
        )
        # updates the correlation parameters (alpha, beta..)
        pars_cor = (config["jacobian"] @ pars_cor.flatten() + config["center"]).reshape(
            (nzs, npar - 1)
        )
        # additive bias parameters (c_1, c_2)
        pars_c1 = np.array(
            [[block[psfsec, "psf_c1_z%d" % (j + 1)]] for j in range(nzs)]
        )
        pars_c2 = np.array(
            [[block[psfsec, "psf_c2_z%d" % (j + 1)]] for j in range(nzs)]
        )
        pars1 = np.hstack([pars_cor, pars_c1])
        pars2 = np.hstack([pars_cor, pars_c2])
        rr = config["r"]
        dxip0 = generate_delta_xip(pars1, pars2, sys1, sys2)
        if nzs == 1:
            dxip0 = np.tile(dxip0, (ncross, 1))
        dxip = Interp1d(rr, dxip0, bounds_error=False)(theta)
    else:
        raise ValueError("model_type: %s does not support" % config["model_type"])

    ib = 0
    for i in range(n_a):
        for j in range(i, n_b):
            bn = "bin_{}_{}".format(j + 1, i + 1)
            # cut the xip
            block["shear_xi_plus", bn] = block["shear_xi_plus", bn][mskT] + dxip[ib]
            ib += 1
    return 0
