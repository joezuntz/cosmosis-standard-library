import numpy as np
from .common import C1_RHOCRIT


def intrinsic_spectra(z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, Omega_m):
    """ 
    The Kirk, Rassat, Host, Bridle (2011) Linear Alignment model.
    Equations 8 and 9.

    The output P_II and P_GI will be specified at the same (k,z) as the linear power inputs
    """

    #extrapolate our linear power out to high redshift
    z0 = np.where(z_lin==0)[0][0]
    nz = len(z_lin)

    # P_II is actually fixed across redshifts
    f = - Omega_m * C1_RHOCRIT

    # intrinsic-intrinsic term
    P_II = np.zeros_like(P_lin)
    for i in xrange(nz):
        P_II[i,:] = f**2 * P_lin[z0,:]


    #Get nonlinear P(k) at the same sampling as linear
    P_nl_resample = np.zeros_like(P_lin)
    for i in xrange(nz):
        log_P_resample = np.interp(np.log(k_lin), np.log(k_nl), np.log(P_nl[i,:]))
        P_nl_resample[i,:] = np.exp(log_P_resample)
    

    growth = np.zeros_like(P_lin)
    ksmall = np.argmin(k_lin)
    for i in xrange(nz):
        growth[i,:] = (P_lin[i,ksmall] / P_lin[z0,ksmall])**0.5

    P_GI = f * P_lin**0.5 * P_nl_resample**0.5 / growth

    return P_II, P_GI
