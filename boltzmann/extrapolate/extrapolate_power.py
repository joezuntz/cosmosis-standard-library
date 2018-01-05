"""
This module takes linear and non-linear P(k) and extrapolates
them linearly in log-space out to a specified high k_max

"""
from builtins import range
from cosmosis.datablock import option_section, names
import numpy as np
from numpy import log, exp


def linear_extend(x, y, xmin, xmax, nmin, nmax, nfit):
    if xmin < x.min():
        xf = x[:nfit]
        yf = y[:nfit]
        p = np.polyfit(xf, yf, 1)
        xnew = np.linspace(xmin, x.min(), nmin, endpoint=False)
        ynew = np.polyval(p, xnew)
        x = np.concatenate((xnew, x))
        y = np.concatenate((ynew, y))
    if xmax > x.max():
        xf = x[-nfit:]
        yf = y[-nfit:]
        p = np.polyfit(xf, yf, 1)
        xnew = np.linspace(x.max(), xmax, nmax, endpoint=True)
        # skip the first point as it is just the xmax
        xnew = xnew[1:]
        ynew = np.polyval(p, xnew)
        x = np.concatenate((x, xnew))
        y = np.concatenate((y, ynew))
    return x, y


def extrapolate_section(block, section, kmin, kmax, nmin, nmax, npoint):
    # load current values
    k = block[section, "k_h"]
    z = block[section, "z"]
    nk = len(k)
    nz = len(z)
    # load other current values
    k, z, P = block.get_grid(section, "k_h", "z", "p_k")
    # extrapolate
    P_out = []
    for i in range(nz):
        Pi = P[:, i]
        logk, logp = linear_extend(log(k), log(Pi), log(
            kmin), log(kmax), nmin, nmax, npoint)
        P_out.append(exp(logp))

    k = exp(logk)
    P_out = np.dstack(P_out).squeeze()

    block.replace_grid(section, "z", z, "k_h", k, "P_k", P_out.T)


def setup(options):
    kmax = options.get_double(option_section, "kmax")
    kmin = options.get_double(option_section, "kmin", default=1e10)
    nmin = options.get_int(option_section, "nmin", default=50)
    npoint = options.get_int(option_section, "npoint", default=3)
    nmax = options.get_int(option_section, "nmax", default=200)
    return {"kmax": kmax, "kmin": kmin, "nmin": nmin, "nmax": nmax, "npoint": npoint}


def execute(block, config):
    kmin = config['kmin']
    kmax = config['kmax']
    nmin = config['nmin']
    nmax = config['nmax']
    npoint = config['npoint']

    # extrapolate non-linear power
    for section in [names.matter_power_nl, names.matter_power_lin]:
        if block.has_section(section):
            extrapolate_section(block, section, kmin, kmax, nmin, nmax, npoint)
    return 0
