#coding: utf-8
from builtins import range
import scipy.interpolate
import numpy as np
"""
This module generates the effective power P_12 
linear alignment KRHB model, which goes into C_ell
(under the Limber approximation) as :

C_ell = \int X_1(chi) X_2(chi) / chi^2 P_12(k=ell/chi,chi) d chi


"""


def compute_c1_baseline():
    C1_M_sun = 5e-14  # h^-2 M_S^-1 Mpc^3
    M_sun = 1.9891e30  # kg
    Mpc_in_m = 3.0857e22  # meters
    C1_SI = C1_M_sun / M_sun * (Mpc_in_m)**3  # h^-2 kg^-1 m^3
    # rho_crit_0 = 3 H^2 / 8 pi G
    G = 6.67384e-11  # m^3 kg^-1 s^-2
    H = 100  #  h km s^-1 Mpc^-1
    H_SI = H * 1000.0 / Mpc_in_m  # h s^-1
    rho_crit_0 = 3 * H_SI**2 / (8 * np.pi * G)  # h^2 kg m^-3
    f = C1_SI * rho_crit_0
    return f


def resample_power(k_new, k_old, P_old):
    "Log-Linearly resample P1[z,k1] onto k2. z values assumed equal"
    nk = len(k_new)
    nz = P_old.shape[0]
    P_new = np.zeros((nz, nk))
    for i in range(nz):
        if (P_old[i] == 0).all():
            P_new[i] = np.zeros(nk)
        else:
            p_i = np.interp(np.log(k_new), np.log(
                k_old), np.log(abs(P_old[i])))
            P_new[i] = np.exp(p_i) * np.sign(P_old[i][0])
    return P_new


# in units of rho_crit0
C1_RHOCRIT = compute_c1_baseline()
# print "C1_RHOCRIT = ", C1_RHOCRIT


def bridle_king(z_nl, k_nl, P_nl, A, Omega_m):
    # extrapolate our linear power out to high redshift
    z0 = np.where(z_nl == 0)[0][0]
    nz = len(z_nl)

    # P_II is actually fixed across redshifts
    # Which forces b_I to vary
    f = - A * Omega_m * C1_RHOCRIT

    # intrinsic-intrinsic term
    P_II = np.zeros_like(P_nl)
    for i in range(nz):
        P_II[i] = f**2 * P_nl[z0]

    growth = np.zeros_like(z_nl)
    ksmall = np.argmin(k_nl)
    P_GI = np.zeros_like(P_nl)
    for i in range(nz):
        growth = (P_nl[i, ksmall] / P_nl[z0, ksmall])**0.5
        P_GI[i] = f * P_nl[i] / growth

    # Finally calculate the intrinsic and stochastic bias terms from the power spectra
    R1 = P_II / P_nl
    b_I = -1.0 * np.sqrt(R1) * np.sign(A)
    r_I = P_GI / P_II * b_I

    return P_II, P_GI, b_I, r_I, k_nl, z_nl


def bridle_king_corrected(z_nl, k_nl, P_nl, A, Omega_m):
    # What was used in CFHTLens and Maccrann et al.
    # extrapolate our linear power out to high redshift
    z0 = np.where(z_nl == 0)[0][0]
    nz = len(z_nl)

    ksmall = np.argmin(k_nl)

    growth = (P_nl[:, ksmall] / P_nl[z0, ksmall])**0.5

    F = - A * C1_RHOCRIT * Omega_m / growth

    # intrinsic-intrinsic term
    P_II = np.zeros_like(P_nl)

    for i in range(nz):
        P_II[i, :] = F[i]**2 * P_nl[i, :]

    P_GI = np.zeros_like(P_nl)
    for i in range(nz):
        P_GI[i] = F[i] * P_nl[i]

    # Finally calculate the intrinsic and stochastic bias terms from the power spectra
    R1 = P_II / P_nl
    b_I = -1.0 * np.sqrt(R1) * np.sign(A)
    r_I = P_GI / P_II * b_I

    return P_II, P_GI, b_I, r_I, k_nl, z_nl


def linear(z_lin, k_lin, P_lin, A, Omega_m):
    # What was used in CFHTLens and Maccrann et al.
    # extrapolate our linear power out to high redshift
    z0 = np.where(z_lin == 0)[0][0]
    nz = len(z_lin)

    ksmall = np.argmin(k_lin)

    growth = (P_lin[:, ksmall] / P_lin[z0, ksmall])**0.5

    F = - A * C1_RHOCRIT * Omega_m / growth

    # intrinsic-intrinsic term
    P_II = np.zeros_like(P_lin)

    for i in range(nz):
        P_II[i, :] = F[i]**2 * P_lin[i, :]

    P_GI = np.zeros_like(P_lin)
    for i in range(nz):
        P_GI[i] = F[i] * P_lin[i]

    # Finally calculate the intrinsic and stochastic bias terms from the power spectra
    R1 = P_II / P_lin
    b_I = np.sqrt(R1) * -1.0 * A / abs(A)
    r_I = P_GI / P_II * b_I

    return P_II, P_GI, b_I, r_I, k_lin, z_lin


def kirk_rassat_host_bridle_power(z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, A, Omega_m):
    """ 
    The Kirk, Rassat, Host, Bridle (2011) Linear Alignment model.
    Equations 8 and 9.

    C1 and rho_m0 must be in consistent units.
    The output P_II and P_GI will be specified at the same (k,z) as the linear power inputs

    The input f = A*C_1*rho_crit0

    """

    # extrapolate our linear power out to high redshift
    z0 = np.where(z_lin == 0)[0][0]
    nz = len(z_lin)

    # P_II is actually fixed across redshifts
    f = - Omega_m * A * C1_RHOCRIT

    # intrinsic-intrinsic term
    P_II = np.zeros_like(P_lin)
    for i in range(nz):
        P_II[i] = f**2 * P_lin[z0]

    # Make sure the k sampling of the two spectra are contiguous
    ikstart = np.argwhere((k_lin > k_nl[0]) & (k_lin < k_nl[-1])).T[0, 0]
    ikend = np.argwhere((k_lin > k_nl[0]) & (k_lin < k_nl[-1])).T[0, -1]

    P_II = P_II[:, ikstart:ikend]
    P_lin = P_lin[:, ikstart:ikend]
    k_lin = k_lin[ikstart:ikend]

    P_nl_resample = resample_power(k_lin, k_nl, P_nl)

    growth = np.zeros_like(P_lin)
    ksmall = np.argmin(k_lin)
    for i in range(nz):
        growth[i] = (P_lin[i, ksmall] / P_lin[z0, ksmall])**0.5

    P_GI = f * P_lin**0.5 * P_nl_resample**0.5 / growth

    # Finally calculate the IA and stochastic bias terms from the power spectra
    P_II_resample = resample_power(k_nl, k_lin, P_II)
    P_GI_resample = resample_power(k_nl, k_lin, P_GI)

    R1 = P_II_resample / P_nl
    b_I = -1.0 * np.sqrt(R1) * np.sign(A)
    r_I = P_GI_resample / P_II_resample * b_I

    return P_II_resample, P_GI_resample, b_I, r_I, k_nl, z_nl
