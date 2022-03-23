# -*- coding: utf-8 -*-
import numpy as np

def compute_c1_baseline():
    C1_M_sun = 5e-14 # h^-2 M_S^-1 Mpc^3
    M_sun = 1.9891e30 # kg
    Mpc_in_m = 3.0857e22 # meters
    C1_SI = C1_M_sun / M_sun * (Mpc_in_m)**3  # h^-2 kg^-1 m^3
    #rho_crit_0 = 3 H^2 / 8 pi G
    G = 6.67384e-11 #m^3 kg^-1 s^-2
    H = 100 # h km s^-1 Mpc^-1
    H_SI = H * 1000.0 / Mpc_in_m  # h s^-1
    rho_crit_0 = 3 * H_SI**2 / (8*np.pi*G)  #  h^2 kg m^-3
    f = C1_SI * rho_crit_0
    return f

C1_RHOCRIT = compute_c1_baseline()

def resample_power(P1, P2, k1, k2):
    "Linearly resample P2 into the k values from P1. z values assumed equal"
    P_resampled = np.zeros_like(P1)
    nz = P1.shape[0]
    for i in xrange(nz):
        p_i = np.interp(np.log(k1), np.log(k2), np.log(abs(P2[i])))
        P_resampled[i] = np.exp(p_i) * np.sign(P2[i,0])
    return P_resampled

