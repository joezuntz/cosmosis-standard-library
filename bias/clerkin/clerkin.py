import numpy as np


def gtd_bias(z, growth, alpha, b0, c):
    b = c + (b0 - c) / growth**alpha
    return b


def q_bias(k, Q, A):
    return (1 + Q * k**2) / (1 + A * k)


def make_grids(k, z):
    K = np.tile(k[:, None], z.size)
    Z = np.tile(z[:, None], k.size).T
    return K, Z


def q_model(k, z, Q, A):
    # Make 2D versions of k,z arrays for convenience
    K, Z = make_grids(k, z)
    bias = q_bias(K, Q, A)
    return bias


def gtd_model(k, z, z_growth, growth, alpha, b0, c):
    K, Z = make_grids(k, z)
    D = np.interp(z, z_growth, growth)
    D = np.tile(D[:, None], k.size).T
    bias = gtd_bias(Z, D, alpha, b0, c)
    return bias


def gtd_q_model(k, z, z_growth, growth, alpha, b0, c, Q, A):
    K, Z = make_grids(k, z)
    bias_k = q_bias(K, Q, A)
    bias_z = gtd_model(k, z, z_growth, growth, alpha, b0, c)
    bias = bias_k * bias_z
    return bias
