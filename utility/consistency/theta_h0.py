import camb
import numpy as np
import scipy.integrate
import scipy.optimize


def H0_to_theta(hubble, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2, omega_k, final=False):
    h = hubble / 100.
    if np.isnan(omnuh2):
        omnuh2 = omega_nu * h**2
    if np.isnan(omega_m):
        omega_m = ommh2 / h**2
    if np.isnan(ombh2):
        ombh2 = omega_b * h**2
    if np.isnan(omega_lambda):
        omega_lambda = omlamh2 / h**2
    if np.isnan(omch2):
        omch2 = omega_c * h**2
    if np.isnan(omega_m):
        omega_m = omega_c + omega_b + omega_nu
    if np.isnan(omega_m):
        omega_m = (omch2 + ombh2 + omnuh2) / h**2
    if np.isnan(omega_lambda):
        omega_lambda = 1 - omega_k - omega_m

    if np.isnan([omnuh2, hubble, omega_m, omega_lambda, ombh2]).any():
        return np.nan
    
    mnu = 93.14 * omnuh2

    original_feedback_level = camb.config.FeedbackLevel

    try:
        camb.set_feedback_level(0)
        p = camb.CAMBparams()
        p.set_cosmology(ombh2 = ombh2,
                        omch2 = omch2,
                        omk = omega_k,
                        mnu = mnu,
                        H0=hubble)
        r = camb.get_background(p)
        theta = r.cosmomc_theta()
    except:
        theta = np.nan
    finally:
        camb.config.FeedbackLevel = original_feedback_level
    
    return theta


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
    omega_k = params.get('omega_k', np.nan)

    return 100 * H0_to_theta(hubble, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2, omega_k)


def theta_to_H0(theta, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2, omega_k):
    get_theta = lambda H0: H0_to_theta(H0, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2, omega_k)
    t = get_theta(70.0)

    if np.isnan(t):
        return np.nan

    H0_min = 20.0
    H0_max = 120.0

    try:
        get_theta(H0_min)
    except ValueError:
        H0_min = 40.0

    try:
        get_theta(H0_max)
    except ValueError:
        H0_max = 100.0

    f = lambda H0: get_theta(H0) - theta
    H0 = scipy.optimize.brentq(f, H0_min, H0_max, rtol=5e-5)

    return H0


def theta_to_H0_interface(params):  
    theta = params['cosmomc_theta'] / 100
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
    omegak = params.get('omega_k', np.nan)

    return theta_to_H0(theta, omega_nu, omnuh2, omega_m, ommh2, omega_c, omch2, omega_b, ombh2, omega_lambda, omlamh2, omegak)
