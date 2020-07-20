from astropy.cosmology import LambdaCDM
from astropy import units, constants
import numpy as np
import scipy.integrate
import scipy.optimize

def dtauda(a, cosmo):
    z = 1.0 / a - 1
    H = cosmo.H(z)
    return  (constants.c / (H * a**2)).to("Mpc").value

def dsound_da(a, cosmo):
    ombh2 = cosmo.Ob0 * cosmo.h**2
    R = 3.0e4 * a * ombh2
    c_s = 1.0 / np.sqrt(3*(1+R))
    return dtauda(a, cosmo) * c_s

def H0_to_theta(hubble, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2):
    h = hubble / 100.
    if np.isnan(omnuh2):
        omnuh2 = omega_nu * h**2
    if np.isnan(omega_m):
        omega_m = ommh2 / h**2
    if np.isnan(omega_b):
        omega_b = ombh2 / h**2
    if np.isnan(omega_lambda):
        omega_lambda = omlamh2 / h**2
    if np.isnan(omega_c):
        omega_c = omch2 / h**2
    if np.isnan(omega_m):
        omega_m = omega_c + omega_b + omega_nu
    if np.isnan(omega_m):
        omega_m = (omch2 + ombh2 + omnuh2) / h**2

    if np.isnan([omnuh2, hubble, omega_m, omega_lambda, omega_b]).any():
        return np.nan
    
    m_nu = 93.14 * omnuh2 * units.eV
    cosmo = LambdaCDM(hubble, omega_m, omega_lambda, m_nu=m_nu/3, Ob0=omega_b, Tcmb0=2.7255, Neff=3.046)

    ombh2 = cosmo.Ob0 * cosmo.h**2
    omdmh2 = cosmo.Om0 * cosmo.h**2

    zstar = 1048*(1+0.00124*ombh2**(-0.738))*(1+ (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))
    astar = 1 / (1+zstar)
    rs = scipy.integrate.romberg(dsound_da, 1e-8, astar, rtol=5e-5, args=(cosmo,), divmax=10) # in Mpc
    DA = cosmo.angular_diameter_distance(zstar).to("Mpc").value / astar
    return rs / DA


def H0_to_theta_interface(params):
    hubble = params['hubble']
    omega_nu = params.get('omega_nu', np.nan)
    omnuh2 = params.get('omnuh2', np.nan)
    omega_m = params.get('omega_m', np.nan)
    ommh2 = params.get('ommh2', np.nan)
    omega_c = params.get('omega_c', np.nan)
    omch2 = params.get('omch2', np.nan)
    omega_b = params.get('omega_b', np.nan)
    ombh2 = params.get('ombh2', np.nan)
    omega_lambda = params.get('omega_lambda', np.nan)
    omlamh2 = params.get('omlamh2', np.nan)

    return H0_to_theta(hubble, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2)


def theta_to_H0(theta, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2):

    get_theta = lambda H0: H0_to_theta(H0, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2)
    t = get_theta(70.0)

    if np.isnan(t):
        return np.nan

    H0_min = 20.0
    H0_max = 120.0

    try:
        theta_min = get_theta(H0_min)
    except ValueError:
        H0_min = 40.0

    try:
        theta_max = get_theta(H0_max)
    except ValueError:
        H0_max = 100.0

    f = lambda H0: get_theta(H0) - theta
    H0 = scipy.optimize.bisect(f, H0_min, H0_max, rtol=5e-5)
    return H0


def theta_to_H0_interface(params):
    theta = params['cosmomc_theta']
    omega_nu = params.get('omega_nu', np.nan)
    omnuh2 = params.get('omnuh2', np.nan)
    omega_m = params.get('omega_m', np.nan)
    ommh2 = params.get('ommh2', np.nan)
    omega_c = params.get('omega_c', np.nan)
    omch2 = params.get('omch2', np.nan)
    omega_b = params.get('omega_b', np.nan)
    ombh2 = params.get('ombh2', np.nan)
    omega_lambda = params.get('omega_lambda', np.nan)
    omlamh2 = params.get('omlamh2', np.nan)

    return theta_to_H0(theta, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2)
