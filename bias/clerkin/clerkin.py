import numpy as np

def gtd_bias(z, growth, alpha, b0, c):
	b = c + (b0-c) / growth**alpha
	return b

def gtd_model(z1, growth_z, k, z, P, alpha, b0, c):
	D = np.interp(z, z1, growth_z)
	b = gtd_bias(z, D, alpha, b0, c)
	P_out = P.copy()
	nk = len(k)
	for i in xrange(nk):
		P_out[i, :] *= b
	return P_out

def q_bias(k, Q, A):
	return (1+Q*k**2) / (1+A*k)

def q_model(k, z, P, Q, A):
	b = q_bias(k, Q, A)
	P_out = P.copy()
	nz = len(z)
	for i in xrange(nz):
		P_out[:, i] *= b
	return P_out

def gtd_q_model(z1, growth_z, k, z, P, alpha, b0, c, Q, A):
	D = np.interp(z, z1, growth_z)
	b_z = gtd_bias(z, D, alpha, b0, c)
	b_k = q_bias(k, Q, A)
	P_out = P.copy()
	nk = len(k)
	nz = len(z)
	for i in xrange(nk):
		P_out[i, :] *= b_z
	for i in xrange(nz):
		P_out[:, i] *= b_k
	return P_out
