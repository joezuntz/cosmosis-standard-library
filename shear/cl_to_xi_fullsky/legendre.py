from __future__ import print_function
from builtins import range
import numpy as np
from scipy.special import lpn

PI=np.pi

def sin_filter(ell_max, ell_right):
    ells = np.arange(ell_max+1)
    y = (ell_max-ells)/(ell_max-ell_right)
    return np.where( ells>ell_right, y - np.sin(2*PI*y)/2/PI, np.ones(len(ells)) )

def get_N_ell(ell):
    """N_ell (eq. 2.16)"""
    ell = np.atleast_1d(ell)
    N_ell = np.sqrt(2. / (ell - 1) / ell / (ell + 1) / (ell + 2))
    N_ell[ell < 2] = 0.
    N_ell[ell == 2] = np.sqrt(2. / (4 * 3 * 2))
    return N_ell

def apply_filter(ell_max, high_l_filter, legfacs):
    f = sin_filter( ell_max, ell_max*high_l_filter)
    return legfacs * np.tile(f, (legfacs.shape[0],1))

def get_F_theta_l(theta_vals, ell_max, cl2xi_type, high_l_filter=0.75):
    ells = np.arange(ell_max+1)
    if cl2xi_type == "00":
        leg_factors = get_legfactors_00(ells, theta_vals)
    elif cl2xi_type == "02+":
        leg_factors = get_legfactors_02(ells, theta_vals)
    elif cl2xi_type in ["22+","22-"]:
        gp, gm = precomp_GpGm(ells, theta_vals)
        N_ell = get_N_ell(ells)
        leg_factors_p, leg_factors_m = get_legfactors_22(ells, theta_vals)
        leg_factors = leg_factors_p if ( cl2xi_type=="22+" ) else leg_factors_m
    else:
        raise ValueError("cl2xi type %s not recognised (should be one of 00, 02+, 22+, 22-)"%cl2xi_type)
    if high_l_filter is not None:
        f = sin_filter( ell_max, ell_max*high_l_filter )
        print("filter:",f)
        leg_factors *= np.tile(f, (len(theta_vals), 1) )

    return leg_factors

def get_legfactors_00(ells, thetas):
    n_ell, n_theta = len(ells), len(thetas)
    legfacs = np.zeros((n_theta, n_ell))
    for it, t in enumerate(thetas):
        legfacs[it] = (2 * ells + 1) * lpn(ells[-1], np.cos(t))[0] / 4. / PI
    legfacs[:, 0] = 0.
    return legfacs

def get_legfactors_02(ells, thetas):
    n_ell, n_theta = len(ells), len(thetas)
    legfacs = np.zeros((n_theta, n_ell))
    ell_factor = np.zeros(len(ells))
    ell_factor[1:] = (2 * ells[1:] + 1) / 4. / PI / ells[1:] / (ells[1:] + 1)
    for it, t in enumerate(thetas):
        P2l = P2l_rec_norm(ells, np.cos(t))
        legfacs[it] = P2l * ell_factor
    return legfacs

def get_legfactors_22(ells, thetas):
    gp, gm = precomp_GpGm(ells, thetas)
    N_ell = get_N_ell(ells)
    prefac = ((2 * ells + 1) / 2. / PI) * N_ell * N_ell / 2
    leg_factors_p = prefac * (gp + gm) 
    leg_factors_m = prefac * (gp - gm) 
    return ( leg_factors_p, leg_factors_m )

def P2l_rec(ells, cost):
    """Calculate P2l using recurrence relation"""
    P22 = 3 * (1 - cost**2)
    P23 = 15 * cost * (1 - cost**2)
    P2l = np.zeros(len(ells))
    P2l[0] = 0.
    P2l[1] = 0.
    P2l[2] = P22
    P2l[3] = P23
    for ell in ells[4:]:
        # print ell, P2l[ell-1], P2l[ell-2]
        P2l[ell] = ((2 * ell - 1) * cost * P2l[ell - 1] -
                    (ell + 2 - 1) * P2l[ell - 2]) / (ell - 2)
    return P2l

def P2l_norm_prefac(ell):
    return np.sqrt((2. * ell + 1) / ((ell - 1.) * ell * (ell + 1.) * (ell + 2.)) / 4. / PI)

def P2l_rec_norm(ells, cost):
    """Calculate P2l using recurrence relation for normalised P2l"""
    P22 = 3. * (1. - cost**2)
    P23 = 15. * cost * (1. - cost**2)
    P2l = np.zeros(len(ells))
    P2l[0] = 0.
    P2l[1] = 0.
    P2l[2] = P22
    P2l[3] = P23
    P2l_norm = np.copy(P2l)
    P2l_norm[2] *= P2l_norm_prefac(2)
    P2l_norm[3] *= P2l_norm_prefac(3)
    for ell in ells[4:]:
        # print ell, P2l[ell-1], P2l[ell-2]
        a = np.sqrt((4 * ell**2 - 1.) / (ell**2 - 4))
        b = cost * P2l_norm[ell - 1]
        c = np.sqrt(((ell - 1.)**2 - 4) /
                    (4 * (ell - 1.)**2 - 1)) * P2l_norm[ell - 2]
        # print a,b,c
        P2l_norm[ell] = a * (b - c)
        # print ell, P2l_norm[ell], P2l_norm_prefac(ell)
        P2l[ell] = P2l_norm[ell] / P2l_norm_prefac(ell)
    return P2l

def precomp_GpGm(ells, thetas):
    """G+/- eqn. 2.26 and 2.27 from astro-ph/9611125v1"""
    n_ell, n_theta = len(ells), len(thetas)
    P_m_l = np.zeros((n_theta, n_ell))
    P_m_lminus1 = np.zeros_like(P_m_l)
    costs = np.cos(thetas)
    sints = np.sin(thetas)
    for it in range(n_theta):
        P_m_l[it] = P2l_rec_norm(ells, costs[it])
        # for il in range(n_ell):
        #    P_m_l[it,il] = sp.lpmv(2, ells[il], costs[it])
    P_m_lminus1[:, 1:] = P_m_l[:, :-1]
    ELLS, THETAS = np.meshgrid(ells, thetas)
    COSTS, SINTS = np.cos(THETAS), np.sin(THETAS)
    G_plus = -((ELLS - 4) / SINTS**2 + 0.5 * ELLS * (ELLS - 1)) * \
        P_m_l + (ELLS + 2) * COSTS * P_m_lminus1 / SINTS**2
    G_minus = 2 * ((ELLS - 1) * COSTS * P_m_l -
                   (ELLS + 2) * P_m_lminus1) / SINTS**2
    G_plus[:, 0] = 0.
    G_minus[:, 0] = 0.
    return G_plus, G_minus

def G_plus_minus_l2(ells, theta, s_switch=100.):
    P_m_l = np.zeros_like(ells)
    P_m_lminus1 = np.zeros_like(P_m_l)
    cost = np.cos(theta)
    sint = np.sin(theta)
    do_full = (ells / theta < s_switch) | (ells < 20)
    for ell in ells[do_full]:
        P_m_l[ell] = sp.lpmv(2, ell, cost)
    P_m_l[~do_full] = (ells[~do_full] - 1) * ells[~do_full] * (ells[~do_full] + 1) * (
        ells[~do_full] + 2) * sp.jn(2, ells[~do_full] * theta) / ells[~do_full]**2
    print('theta = %f, switching to bessel function for %d<ell<%d' % (np.degrees(theta) * 60., ells[~do_full][0], ells[~do_full][-1]))
    P_m_lminus1[1:] = P_m_l[:-1]
    G_plus = -((ells - 4) / sint**2 + 0.5 * ells * (ells - 1)) * \
        P_m_l + (ells + 2) * cost * P_m_lminus1 / sint**2
    G_minus = 2 * ((ells - 1) * cost * P_m_l -
                   (ells + 2) * P_m_lminus1) / sint**2
    G_plus[0] = 0.
    G_minus[0] = 0.
    return G_plus, G_minus
