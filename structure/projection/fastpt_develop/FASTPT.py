"""
	FASTPT is a numerical algorithm to calculate
	1-loop contributions to the matter power spectrum
	and other integrals of a similar type.
	The method is presented in papers arXiv:1603.04826 and arXiv:1609.05978
	Please cite these papers if you are using FASTPT in your research.

	Joseph E. McEwen (c) 2016
	mcewen.24@osu.edu

	Xiao Fang
	fang.307@osu.edu

	Jonathan A. Blazek
	blazek.35@osu.edu


	FFFFFFFF    A           SSSSSSSSS   TTTTTTTTTTTTTT             PPPPPPPPP    TTTTTTTTTTTT
	FF     	   A A         SS                 TT                   PP      PP        TT
	FF        A   A        SS                 TT                   PP      PP        TT
	FFFFF    AAAAAAA        SSSSSSSS          TT       ==========  PPPPPPPPP         TT
	FF      AA     AA              SS         TT                   PP                TT
	FF     AA       AA             SS         TT                   PP                TT
	FF    AA         AA    SSSSSSSSS          TT                   PP                TT


	The FASTPT class is the workhorse of the FASTPT algorithm.
	This class calculates integrals of the form:

	\int \frac{d^3q}{(2 \pi)^3} K(q,k-q) P(q) P(|k-q|)

	\int \frac{d^3q_1}{(2 \pi)^3} K(\hat{q_1} \dot \hat{q_2},\hat{q_1} \dot \hat{k}, \hat{q_2} \dot \hat{k}, q_1, q_2) P(q_1) P(|k-q_1|)

"""
from __future__ import division
from __future__ import print_function

from info import __version__

import numpy as np
from numpy.fft import fft, ifft, rfft, irfft, fftfreq
from numpy import exp, log, log10, cos, sin, pi, cosh, sinh, sqrt
from scipy.special import gamma
from scipy.signal import fftconvolve
import scipy.integrate as integrate
from fastpt_extr import p_window, c_window, pad_left, pad_right
from matter_power_spt import P_13_reg, Y1_reg_NL, Y2_reg_NL
from initialize_params import scalar_stuff, tensor_stuff
from IA_tt import IA_tt
from IA_ABD import IA_A, IA_DEE, IA_DBB, P_IA_B
from IA_ta import IA_deltaE1, P_IA_deltaE2, IA_0E0E, IA_0B0B
from OV import OV
from kPol import kPol
from RSD import RSDA, RSDB
import RSD_ItypeII
from P_extend import k_extend
import FASTPT_simple as fastpt_simple

## WHEN DOES THE IMPORT STEP OCCUR? DO WE WANT TO MOVE SOME OF THESE TO THE INITIALIZATION BLOCK TO SAVE TIME ON LIGHT RUNS?


log2 = log(2.)


class FASTPT:

    def __init__(self, k, nu=None, to_do=None, param_mat=None, low_extrap=None, high_extrap=None, n_pad=None,
                 verbose=False):

        """ inputs:
				* k grid
				* the to_do list: e.g. one_loop density density , bias terms, ...
				* low_extrap is the call to extrapolate the power spectrum to lower k-values,
					this helps with edge effects
				* n_pad is the number of zeros to add to both ends of the array. This helps with
					edge effects.
				* verbose is to turn on verbose settings.
		"""

        # if no to_do list is given, default to fastpt_simple SPT case
        if to_do is None:
            if verbose:
                print(
                    'Note: You are using an earlier call structure for FASTPT. Your code will still run correctly, calling FASTPT_simple. See user manual.')
            if nu is None:  # give a warning if nu=None that a default value is being used.
                print('WARNING: No value for nu is given. FASTPT_simple is being called with a default of nu=-2')
                nu = -2  # this is the default value for P22+P13 and bias calculation
            self.pt_simple = fastpt_simple.FASTPT(k, nu, param_mat=param_mat, low_extrap=low_extrap,
                                                  high_extrap=high_extrap, n_pad=n_pad, verbose=verbose)
            return None
        # Exit initialization here, since fastpt_simple performs the various checks on the k grid and does extrapolation.

        # check for log spacing
        print('Initializing k-grid quantities...')
        dk = np.diff(np.log(k))
        # dk_test=np.ones_like(dk)*dk[0]
        delta_L = (log(k[-1]) - log(k[0])) / (k.size - 1)
        dk_test = np.ones_like(dk) * delta_L

        log_sample_test = 'ERROR! FASTPT will not work if your in put (k,Pk) values are not sampled evenly in log space!'
        np.testing.assert_array_almost_equal(dk, dk_test, decimal=4, err_msg=log_sample_test, verbose=False)

        if verbose:
            print('the minumum and maximum inputed log10(k) are :', np.min(np.log10(k)), np.max(np.log10(k)))
            print('the grid spacing Delta log (k) is', (log(np.max(k)) - log(np.min(k))) / (k.size - 1))
            print('number of input k points are', k.size)
            print('the power spectrum is extraplated to log10(k_min)=', low_extrap)
            print('the power spectrum is extraplated to log10(k_max)=', high_extrap)
            print('the power spectrum has ', n_pad, ' zeros added to both ends of the power spectrum')

        self.k_original = k
        self.extrap = False
        if low_extrap is not None or high_extrap is not None:
            self.EK = k_extend(k, low_extrap, high_extrap)
            k = self.EK.extrap_k()
            self.extrap = True

        self.low_extrap = low_extrap
        self.high_extrap = high_extrap

        self.k_old = k

        # print(self.k_old.size, 'k size')
        # size of input array must be an even number
        if k.size % 2 != 0:
            raise ValueError('Input array must contain an even number of elements.')
        # can we just force the extrapolation to add an element if we need one more? how do we prevent the extrapolation from giving us an odd number of elements? is that hard coded into extrap? or just trim the lowest k value if there is an odd numebr and no extrapolation is requested.

        if n_pad != None:

            self.id_pad = np.arange(k.size) + n_pad
            d_logk = delta_L
            k_pad = np.log(k[0]) - np.arange(1, n_pad + 1) * d_logk
            k_pad = np.exp(k_pad)
            k_left = k_pad[::-1]

            k_pad = np.log(k[-1]) + np.arange(1, n_pad + 1) * d_logk
            k_right = np.exp(k_pad)
            k = np.hstack((k_left, k, k_right))
            n_pad_check = int(np.log(2) / delta_L) + 1
            if n_pad < n_pad_check:
                print('*** Warning ***')
                print('You should consider increasing your zero padding to at least ', n_pad_check)
                print('to ensure that the minimum k_output is > 2k_min in the FASTPT universe.')
                print('k_min in the FASTPT universe is ', k[0], ' while k_min_input is ', self.k_old[0])

        self.k = k
        self.k_size = k.size
        # self.scalar_nu=-2
        self.N = k.size

        # define eta_m and eta_n=eta_m
        omega = 2 * pi / (float(self.N) * delta_L)
        self.m = np.arange(-self.N // 2, self.N // 2 + 1)
        self.eta_m = omega * self.m

        self.verbose = verbose
        self.n_pad = n_pad

        # define l and tau_l
        self.n_l = self.m.size + self.m.size - 1
        self.l = np.arange(-self.n_l // 2 + 1, self.n_l // 2 + 1)
        self.tau_l = omega * self.l

        self.dd_do = False
        self.dd_bias_do = False
        self.IA_tt_do = False
        self.IA_ta_do = False
        self.IA_mix_do = False
        self.OV_do = False
        self.kPol_do = False
        self.RSD_do = False

        for entry in to_do:  # convert to_do list to instructions for FAST-PT initialization
            if entry == 'one_loop_dd':
                self.dd_do = True
                continue
            elif entry == 'dd_bias':
                self.dd_do = True
                self.dd_bias_do = True
                continue
            elif entry == 'IA_all' or entry == 'IA':
                self.IA_tt_do = True
                self.IA_ta_do = True
                self.IA_mix_do = True
                continue
            elif entry == 'IA_tt':
                self.IA_tt_do = True
                continue
            elif entry == 'IA_ta':
                self.IA_ta_do = True
                continue
            elif entry == 'IA_mix':
                self.IA_mix_do = True
                continue
            elif entry == 'OV':
                self.OV_do = True
                continue
            elif entry == 'kPol':
                self.kPol_do = True
                continue
            elif entry == 'RSD':
                self.RSD_do = True
                continue
            elif entry == 'IRres':
                self.dd_do = True
                continue
            elif entry == 'all' or entry == 'everything':
                self.dd_do = True
                self.dd_bias_do = True
                self.IA_tt_do = True
                self.IA_ta_do = True
                self.IA_mix_do = True
                self.OV_do = True
                self.kPol_do = True
                self.RSD_do = True
                continue
            else:
                raise ValueError('FAST-PT does not recognize "' + entry + '" in the to_do list.')

        ### INITIALIZATION of k-grid quantities ###
        if self.dd_do:
            nu = -2
            # parameter matrix for 1-loop calculations
            p_mat = np.array([[0, 0, 0, 0], [0, 0, 2, 0], [0, 0, 4, 0], [2, -2, 2, 0],
                              [1, -1, 1, 0], [1, -1, 3, 0], [2, -2, 0, 1]])

            p_mat_lpt = np.array([[0, 0, 0, 0], [0, 0, 2, 0], [2, -2, 2, 0],
                                  [1, -1, 1, 0], [1, -1, 3, 0], [0, 0, 4, 0], [2, -2, 0, 1]])

            self.X_spt = scalar_stuff(p_mat, nu, self.N, self.m, self.eta_m, self.l, self.tau_l)
            self.X_lpt = scalar_stuff(p_mat_lpt, nu, self.N, self.m, self.eta_m, self.l, self.tau_l)

        if self.IA_tt_do:
            hE_tab, hB_tab = IA_tt()
            p_mat_E = hE_tab[:, [0, 1, 5, 6, 7, 8, 9]]
            p_mat_B = hB_tab[:, [0, 1, 5, 6, 7, 8, 9]]

            self.X_IA_E = tensor_stuff(p_mat_E, self.N, self.m, self.eta_m, self.l, self.tau_l)
            self.X_IA_B = tensor_stuff(p_mat_B, self.N, self.m, self.eta_m, self.l, self.tau_l)

        if self.IA_mix_do:
            IA_A_tab = IA_A()
            IA_DEE_tab = IA_DEE()
            IA_DBB_tab = IA_DBB()
            p_mat_A = IA_A_tab[:, [0, 1, 5, 6, 7, 8, 9]]
            p_mat_DEE = IA_DEE_tab[:, [0, 1, 5, 6, 7, 8, 9]]
            p_mat_DBB = IA_DBB_tab[:, [0, 1, 5, 6, 7, 8, 9]]

            self.X_IA_A = tensor_stuff(p_mat_A, self.N, self.m, self.eta_m, self.l, self.tau_l)
            self.X_IA_DEE = tensor_stuff(p_mat_DEE, self.N, self.m, self.eta_m, self.l, self.tau_l)
            self.X_IA_DBB = tensor_stuff(p_mat_DBB, self.N, self.m, self.eta_m, self.l, self.tau_l)

        if self.IA_ta_do:
            IA_deltaE1_tab = IA_deltaE1()
            IA_0E0E_tab = IA_0E0E()
            IA_0B0B_tab = IA_0B0B()
            p_mat_deltaE1 = IA_deltaE1_tab[:, [0, 1, 5, 6, 7, 8, 9]]
            p_mat_0E0E = IA_0E0E_tab[:, [0, 1, 5, 6, 7, 8, 9]]
            p_mat_0B0B = IA_0B0B_tab[:, [0, 1, 5, 6, 7, 8, 9]]
            self.X_IA_deltaE1 = tensor_stuff(p_mat_deltaE1, self.N, self.m, self.eta_m, self.l, self.tau_l)
            self.X_IA_0E0E = tensor_stuff(p_mat_0E0E, self.N, self.m, self.eta_m, self.l, self.tau_l)
            self.X_IA_0B0B = tensor_stuff(p_mat_0B0B, self.N, self.m, self.eta_m, self.l, self.tau_l)

        if self.OV_do:
            # For OV, we can use two different values for
            # nu1=0 and nu2=-2

            OV_tab = OV()
            p_mat = OV_tab[:, [0, 1, 5, 6, 7, 8, 9]]

            self.X_OV = tensor_stuff(p_mat, self.N, self.m, self.eta_m, self.l, self.tau_l)

        if self.kPol_do:
            tab1, tab2, tab3 = kPol()
            p_mat = tab1[:, [0, 1, 5, 6, 7, 8, 9]]
            self.X_kP1 = tensor_stuff(p_mat, self.N, self.m, self.eta_m, self.l, self.tau_l)

            p_mat = tab2[:, [0, 1, 5, 6, 7, 8, 9]]
            self.X_kP2 = tensor_stuff(p_mat, self.N, self.m, self.eta_m, self.l, self.tau_l)

            p_mat = tab3[:, [0, 1, 5, 6, 7, 8, 9]]
            self.X_kP3 = tensor_stuff(p_mat, self.N, self.m, self.eta_m, self.l, self.tau_l)

        if self.RSD_do:
            tabA, self.A_coeff = RSDA()
            p_mat = tabA[:, [0, 1, 5, 6, 7, 8, 9]]
            self.X_RSDA = tensor_stuff(p_mat, self.N, self.m, self.eta_m, self.l, self.tau_l)

            tabB, self.B_coeff = RSDB()
            p_mat = tabB[:, [0, 1, 5, 6, 7, 8, 9]]
            self.X_RSDB = tensor_stuff(p_mat, self.N, self.m, self.eta_m, self.l, self.tau_l)

    ### Top-level functions to output final quantities ###
    def one_loop_dd(self, P, P_window=None, C_window=None):
        nu = -2

        # routine for one-loop spt calculations

        # coefficents for one_loop calculation
        one_loop_coef = np.array(
            [2 * 1219 / 1470., 2 * 671 / 1029., 2 * 32 / 1715., 2 * 1 / 3., 2 * 62 / 35., 2 * 8 / 35., 1 / 3.])

        # get the roundtrip Fourier power spectrum, i.e. P=IFFT[FFT[P]]
        # get the matrix for each J_k component
        Ps, mat = self.J_k_scalar(P, self.X_spt, nu, P_window=P_window, C_window=C_window)

        P22_mat = np.multiply(one_loop_coef, np.transpose(mat))
        P22 = np.sum(P22_mat, 1)
        P13 = P_13_reg(self.k_old, Ps)
        P_1loop = P22 + P13

        if self.dd_bias_do:
            # if dd_bias is in to_do, this function acts like one_loop_dd_bias

            # Quadraric bias Legendre components
            # See eg section B of Baldauf+ 2012 (arxiv: 1201.4827)
            # Note pre-factor convention is not standardized
            # Returns relevant correlations (including contraction factors),
            # but WITHOUT bias values and other pre-factors.
            # Uses standard "full initialization" of J terms
            sig4 = np.trapz(self.k_old ** 3 * Ps ** 2, x=np.log(self.k_old)) / (2. * pi ** 2)
            # sig4 much more accurate when calculated in logk, especially for low-res input.

            Pd1d2 = 2. * (17. / 21 * mat[0, :] + mat[4, :] + 4. / 21 * mat[1, :])
            Pd2d2 = 2. * (mat[0, :])
            Pd1s2 = 2. * (8. / 315 * mat[0, :] + 4. / 15 * mat[4, :] + 254. / 441 * mat[1, :] + 2. / 5 * mat[5,
                                                                                                         :] + 16. / 245 * mat[
                                                                                                                          2,
                                                                                                                          :])
            Pd2s2 = 2. * (2. / 3 * mat[1, :])
            Ps2s2 = 2. * (4. / 45 * mat[0, :] + 8. / 63 * mat[1, :] + 8. / 35 * mat[2, :])
            if self.extrap:
                _, Ps = self.EK.PK_original(Ps)
                _, P_1loop = self.EK.PK_original(P_1loop)
                _, Pd1d2 = self.EK.PK_original(Pd1d2)
                _, Pd2d2 = self.EK.PK_original(Pd2d2)
                _, Pd1s2 = self.EK.PK_original(Pd1s2)
                _, Pd2s2 = self.EK.PK_original(Pd2s2)
                _, Ps2s2 = self.EK.PK_original(Ps2s2)

            return P_1loop, Ps, Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig4

        if self.extrap:
            _, Ps = self.EK.PK_original(Ps)
            _, P_1loop = self.EK.PK_original(P_1loop)

        return P_1loop, Ps

    def one_loop_dd_bias(self, P, P_window=None, C_window=None):
        nu = -2

        # routine for one-loop spt calculations

        # coefficents for one_loop calculation
        one_loop_coef = np.array(
            [2 * 1219 / 1470., 2 * 671 / 1029., 2 * 32 / 1715., 2 * 1 / 3., 2 * 62 / 35., 2 * 8 / 35., 1 / 3.])

        # get the roundtrip Fourier power spectrum, i.e. P=IFFT[FFT[P]]
        # get the matrix for each J_k component
        Ps, mat = self.J_k_scalar(P, self.X_spt, nu, P_window=P_window, C_window=C_window)

        P22_mat = np.multiply(one_loop_coef, np.transpose(mat))
        P22 = np.sum(P22_mat, 1)
        P13 = P_13_reg(self.k_old, Ps)
        P_1loop = P22 + P13

        # Quadraric bias Legendre components
        # See eg section B of Baldauf+ 2012 (arxiv: 1201.4827)
        # Note pre-factor convention is not standardized
        # Returns relevant correlations (including contraction factors),
        # but WITHOUT bias values and other pre-factors.
        # Uses standard "full initialization" of J terms
        sig4 = np.trapz(self.k_old ** 3 * Ps ** 2, x=np.log(self.k_old)) / (2. * pi ** 2)
        Pd1d2 = 2. * (17. / 21 * mat[0, :] + mat[4, :] + 4. / 21 * mat[1, :])
        Pd2d2 = 2. * (mat[0, :])
        Pd1s2 = 2. * (8. / 315 * mat[0, :] + 4. / 15 * mat[4, :] + 254. / 441 * mat[1, :] + 2. / 5 * mat[5,
                                                                                                     :] + 16. / 245 * mat[
                                                                                                                      2,
                                                                                                                      :])
        Pd2s2 = 2. * (2. / 3 * mat[1, :])
        Ps2s2 = 2. * (4. / 45 * mat[0, :] + 8. / 63 * mat[1, :] + 8. / 35 * mat[2, :])

        if self.extrap:
            _, Ps = self.EK.PK_original(Ps)
            _, P_1loop = self.EK.PK_original(P_1loop)
            _, Pd1d2 = self.EK.PK_original(Pd1d2)
            _, Pd2d2 = self.EK.PK_original(Pd2d2)
            _, Pd1s2 = self.EK.PK_original(Pd1s2)
            _, Pd2s2 = self.EK.PK_original(Pd2s2)
            _, Ps2s2 = self.EK.PK_original(Ps2s2)

        #			return P_1loop, Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig4, Ps #original
        return P_1loop, Ps, Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig4  # new,for consistency

    def one_loop_dd_bias_b3nl(self, P, P_window=None, C_window=None):
        nu = -2

        # routine for one-loop spt calculations

        # coefficents for one_loop calculation
        one_loop_coef = np.array(
            [2 * 1219 / 1470., 2 * 671 / 1029., 2 * 32 / 1715., 2 * 1 / 3., 2 * 62 / 35., 2 * 8 / 35., 1 / 3.])

        # get the roundtrip Fourier power spectrum, i.e. P=IFFT[FFT[P]]
        # get the matrix for each J_k component
        Ps, mat = self.J_k_scalar(P, self.X_spt, nu, P_window=P_window, C_window=C_window)

        P22_mat = np.multiply(one_loop_coef, np.transpose(mat))
        P22 = np.sum(P22_mat, 1)
        P13 = P_13_reg(self.k_old, Ps)
        P_1loop = P22 + P13

        sig4 = np.trapz(self.k_old ** 3 * Ps ** 2, x=np.log(self.k_old)) / (2. * pi ** 2)
        Pd1d2 = 2. * (17. / 21 * mat[0, :] + mat[4, :] + 4. / 21 * mat[1, :])
        Pd2d2 = 2. * (mat[0, :])
        Pd1s2 = 2. * (8. / 315 * mat[0, :] + 4. / 15 * mat[4, :] + 254. / 441 * mat[1, :] + 2. / 5 * mat[5,
                                                                                                     :] + 16. / 245 * mat[
                                                                                                                      2,
                                                                                                                      :])
        Pd2s2 = 2. * (2. / 3 * mat[1, :])
        Ps2s2 = 2. * (4. / 45 * mat[0, :] + 8. / 63 * mat[1, :] + 8. / 35 * mat[2, :])
        sig3nl = Y1_reg_NL(self.k_old, Ps)

        if self.extrap:
            _, Ps = self.EK.PK_original(Ps)
            _, P_1loop = self.EK.PK_original(P_1loop)
            _, Pd1d2 = self.EK.PK_original(Pd1d2)
            _, Pd2d2 = self.EK.PK_original(Pd2d2)
            _, Pd1s2 = self.EK.PK_original(Pd1s2)
            _, Pd2s2 = self.EK.PK_original(Pd2s2)
            _, Ps2s2 = self.EK.PK_original(Ps2s2)
            _, sig3nl = self.EK.PK_original(sig3nl)

        #			return P_1loop, Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig4, Ps #original
        return P_1loop, Ps, Pd1d2, Pd2d2, Pd1s2, Pd2s2, Ps2s2, sig3nl, sig4  # new,for consistency

    def one_loop_dd_bias_lpt_NL(self, P, P_window=None, C_window=None):
        nu_arr = -2

        # get the roundtrip Fourier power spectrum, i.e. P=IFFT[FFT[P]]
        # get the matrix for each J_k component
        Ps, mat = self.J_k_scalar(P, self.X_lpt, nu_arr, P_window=P_window, C_window=C_window)

        [j000, j002, j2n22, j1n11, j1n13, j004, j2n20] = [mat[0, :], mat[1, :], mat[2, :], mat[3, :], mat[4, :],
                                                          mat[5, :], mat[6, :]]

        P22 = 2. * ((1219. / 1470.) * j000 + (671. / 1029.) * j002 + (32. / 1715.) * j004 + (1. / 3.) * j2n22 + (
                62. / 35.) * j1n11 + (8. / 35.) * j1n13 + (1. / 6.) * j2n20)

        sig4 = np.trapz(self.k_old ** 3 * Ps ** 2, x=np.log(self.k_old)) / (2. * pi ** 2)

        X1 = ((144. / 245.) * j000 - (176. / 343.) * j002 - (128. / 1715.) * j004 + (16. / 35.) * j1n11 - (
                16. / 35.) * j1n13)
        X2 = ((16. / 21.) * j000 - (16. / 21.) * j002 + (16. / 35.) * j1n11 - (16. / 35.) * j1n13)
        X3 = (50. / 21.) * j000 + 2. * j1n11 - (8. / 21.) * j002
        X4 = (34. / 21.) * j000 + 2. * j1n11 + (8. / 21.) * j002
        X5 = j000

        Y1 = Y1_reg_NL(self.k_old, Ps)
        Y2 = Y2_reg_NL(self.k_old, Ps)

        Pb1L = X1 + Y1
        Pb1L_2 = X2 + Y2
        Pb1L_b2L = X3
        Pb2L = X4
        Pb2L_2 = X5

        if self.extrap:
            _, Ps = self.EK.PK_original(Ps)
            # _, P_1loop=self.EK.PK_original(P_1loop)

            _, Pb1L = self.EK.PK_original(Pb1L)
            _, Pb1L_2 = self.EK.PK_original(Pb1L_2)
            _, Pb1L_b2L = self.EK.PK_original(Pb1L_b2L)
            _, Pb2L = self.EK.PK_original(Pb2L)
            _, Pb2L_2 = self.EK.PK_original(Pb2L_2)
            _, X1 = self.EK.PK_original(X1)
            _, X2 = self.EK.PK_original(X2)
            _, X3 = self.EK.PK_original(X3)
            _, X4 = self.EK.PK_original(X4)
            _, X5 = self.EK.PK_original(X5)
            _, Y1 = self.EK.PK_original(Y1)
            _, Y2 = self.EK.PK_original(Y2)

        return Ps, Pb1L, Pb1L_2, Pb1L_b2L, Pb2L, Pb2L_2, sig4

    def IA_tt(self, P, P_window=None, C_window=None):

        P_E, A = self.J_k_tensor(P, self.X_IA_E, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_E = self.EK.PK_original(P_E)

        P_B, A = self.J_k_tensor(P, self.X_IA_B, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_B = self.EK.PK_original(P_B)
        return 2. * P_E, 2. * P_B

    ## eq 21 EE; eq 21 BB

    def IA_mix(self, P, P_window=None, C_window=None):

        P_A, A = self.J_k_tensor(P, self.X_IA_A, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_A = self.EK.PK_original(P_A)

        P_Btype2 = P_IA_B(self.k_original, P)

        P_DEE, A = self.J_k_tensor(P, self.X_IA_DEE, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_DEE = self.EK.PK_original(P_DEE)

        P_DBB, A = self.J_k_tensor(P, self.X_IA_DBB, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_DBB = self.EK.PK_original(P_DBB)

        return 2 * P_A, P_Btype2, 2 * P_DEE, 2 * P_DBB

    ## eq 18; eq 19; eq 27 EE; eq 27 BB

    def IA_ta(self, P, P_window=None, C_window=None):

        P_deltaE1, A = self.J_k_tensor(P, self.X_IA_deltaE1, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_deltaE1 = self.EK.PK_original(P_deltaE1)

        P_deltaE2 = P_IA_deltaE2(self.k_original, P)

        P_0E0E, A = self.J_k_tensor(P, self.X_IA_0E0E, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_0E0E = self.EK.PK_original(P_0E0E)

        P_0B0B, A = self.J_k_tensor(P, self.X_IA_0B0B, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P_0B0B = self.EK.PK_original(P_0B0B)

        return 2. * P_deltaE1, 2. * P_deltaE2, P_0E0E, P_0B0B

    ## eq 12 (line 2); eq 12 (line 3); eq 15 EE; eq 15 BB

    def OV(self, P, P_window=None, C_window=None):
        P, A = self.J_k_tensor(P, self.X_OV, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P = self.EK.PK_original(P)
        return P * (2 * pi) ** 2

    def kPol(self, P, P_window=None, C_window=None):
        P1, A = self.J_k_tensor(P, self.X_kP1, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P1 = self.EK.PK_original(P1)

        P2, A = self.J_k_tensor(P, self.X_kP2, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P2 = self.EK.PK_original(P2)

        P3, A = self.J_k_tensor(P, self.X_kP3, P_window=P_window, C_window=C_window)
        if self.extrap:
            _, P3 = self.EK.PK_original(P3)
        return P1 / (80 * pi ** 2), P2 / (160 * pi ** 2), P3 / (80 * pi ** 2)

    def RSD_components(self, P, f, P_window=None, C_window=None):

        _, A = self.J_k_tensor(P, self.X_RSDA, P_window=P_window, C_window=C_window)

        A1 = np.dot(self.A_coeff[:, 0], A) + f * np.dot(self.A_coeff[:, 1], A) + f ** 2 * np.dot(self.A_coeff[:, 2], A)
        A3 = np.dot(self.A_coeff[:, 3], A) + f * np.dot(self.A_coeff[:, 4], A) + f ** 2 * np.dot(self.A_coeff[:, 5], A)
        A5 = np.dot(self.A_coeff[:, 6], A) + f * np.dot(self.A_coeff[:, 7], A) + f ** 2 * np.dot(self.A_coeff[:, 8], A)

        _, B = self.J_k_tensor(P, self.X_RSDB, P_window=P_window, C_window=C_window)

        B0 = np.dot(self.B_coeff[:, 0], B) + f * np.dot(self.B_coeff[:, 1], B) + f ** 2 * np.dot(self.B_coeff[:, 2], B)
        B2 = np.dot(self.B_coeff[:, 3], B) + f * np.dot(self.B_coeff[:, 4], B) + f ** 2 * np.dot(self.B_coeff[:, 5], B)
        B4 = np.dot(self.B_coeff[:, 6], B) + f * np.dot(self.B_coeff[:, 7], B) + f ** 2 * np.dot(self.B_coeff[:, 8], B)
        B6 = np.dot(self.B_coeff[:, 9], B) + f * np.dot(self.B_coeff[:, 10], B) + f ** 2 * np.dot(self.B_coeff[:, 11],
                                                                                                  B)

        if self.extrap:
            _, A1 = self.EK.PK_original(A1)
            _, A3 = self.EK.PK_original(A3)
            _, A5 = self.EK.PK_original(A5)
            _, B0 = self.EK.PK_original(B0)
            _, B2 = self.EK.PK_original(B2)
            _, B4 = self.EK.PK_original(B4)
            _, B6 = self.EK.PK_original(B6)

        P_Ap1 = RSD_ItypeII.P_Ap1(self.k_original, P, f)
        P_Ap3 = RSD_ItypeII.P_Ap3(self.k_original, P, f)
        P_Ap5 = RSD_ItypeII.P_Ap5(self.k_original, P, f)

        return A1, A3, A5, B0, B2, B4, B6, P_Ap1, P_Ap3, P_Ap5

    def RSD_ABsum_components(self, P, f, P_window=None, C_window=None):

        A1, A3, A5, B0, B2, B4, B6, P_Ap1, P_Ap3, P_Ap5 = self.RSD_components(P, f, P_window, C_window)
        ABsum_mu2 = self.k_original * f * (A1 + P_Ap1) + (f * self.k_original) ** 2 * B0
        ABsum_mu4 = self.k_original * f * (A3 + P_Ap3) + (f * self.k_original) ** 2 * B2
        ABsum_mu6 = self.k_original * f * (A5 + P_Ap5) + (f * self.k_original) ** 2 * B4
        ABsum_mu8 = (f * self.k_original) ** 2 * B6

        return ABsum_mu2, ABsum_mu4, ABsum_mu6, ABsum_mu8

    def RSD_ABsum_mu(self, P, f, mu_n, P_window=None, C_window=None):
        ABsum_mu2, ABsum_mu4, ABsum_mu6, ABsum_mu8 = self.RSD_ABsum_components(P, f, P_window, C_window)
        ABsum = ABsum_mu2 * mu_n ** 2 + ABsum_mu4 * mu_n ** 4 + ABsum_mu6 * mu_n ** 6 + ABsum_mu8 * mu_n ** 8
        return ABsum

    def IRres(self, P, L=0.2, h=0.67, rsdrag=135, P_window=None, C_window=None):
        # based on script by M. Ivanov. See arxiv:1605.02149, eq 7.4

        # put this function in the typical fast-pt format, with minimal additional function calls.
        # We can also include a script to do additional things: e.g read in r_BAO from class/camb output
        # or calculate r_BAO from cosmological params.
        from scipy import interpolate
        k = self.k_original
        rbao = h * rsdrag * 1.027  # linear BAO scale
        # set up splining to create P_nowiggle
        kmin = k[0]
        kmax = k[-1]
        knode1 = pi / rbao
        knode2 = 2 * pi / rbao
        klogleft = np.arange(log(kmin), log(3.e-3), 0.2)
        klogright = np.arange(log(0.6), log(kmax), 0.085)
        klogright = np.hstack((log(knode1), log(knode2), klogright))
        kloglist = np.concatenate((klogleft, klogright))
        klist = np.exp(kloglist)

        # how to deal with extended k and P? which values should be used here? probably the extended versions?
        plin = interpolate.InterpolatedUnivariateSpline(k, P)
        logPs = np.log(plin(klist))
        logpsmooth = interpolate.InterpolatedUnivariateSpline(kloglist, logPs)

        def psmooth(x):
            return exp(logpsmooth(log(x)))

        def pw(x):
            return plin(x) - psmooth(x)

        # compute Sigma^2 and the tree-level IR-resummed PS

        Sigma = integrate.quad(lambda x: (4 * pi) * psmooth(x) * (
                1 - 3 * (2 * rbao * x * cos(x * rbao) + (-2 + rbao ** 2 * x ** 2) * sin(rbao * x)) / (
                x * rbao) ** 3) / (3 * (2 * pi) ** 3), kmin, L)[0]

        # speed up by using trap rule integration?
        # change to integration over log-k(?):
        # 		Sigma = integrate.quad(lambda x: x*(4*pi)*psmooth(x)*(1-3*(2*rbao*x*cos(x*rbao)+(-2+rbao**2*x**2)*sin(rbao*x))/(x*rbao)**3)/(3*(2*pi)**3), np.log(kmin), np.log(L))[0]
        def presum(x):
            return psmooth(x) + pw(x) * exp(-x ** 2 * Sigma)

        P_in = presum(k)
        out_1loop = self.one_loop_dd(P_in, P_window=P_window, C_window=C_window)[0]
        # p1loop = interpolate.InterpolatedUnivariateSpline(k,out_1loop) # is this necessary? out_1loop should already be at the correct k spacing
        return psmooth(k) + out_1loop + pw(k) * exp(-k ** 2 * Sigma) * (1 + Sigma * k ** 2)

    ######################################################################################
    ### functions that use the older version structures. ###
    def one_loop(self, P, P_window=None, C_window=None):

        return self.pt_simple.one_loop(P, P_window=P_window, C_window=C_window)

    def P_bias(self, P, P_window=None, C_window=None):

        return self.pt_simple.P_bias(P, P_window=P_window, C_window=C_window)

    ######################################################################################
    ### Core functions used by top-level functions ###
    def J_k_scalar(self, P_in, X, nu, P_window=None, C_window=None):

        pf, p, g_m, g_n, two_part_l, h_l = X

        if self.low_extrap is not None:
            P_in = self.EK.extrap_P_low(P_in)

        if self.high_extrap is not None:
            P_in = self.EK.extrap_P_high(P_in)

        P_b = P_in * self.k_old ** (-nu)

        if self.n_pad is not None:
            P_b = np.pad(P_b, pad_width=(self.n_pad, self.n_pad), mode='constant', constant_values=0)

        c_m_positive = rfft(P_b)
        # We always filter the Fourier coefficients, so the last element is zero.
        # But in case someone does not filter, divide the end point by two
        c_m_positive[-1] = c_m_positive[-1] / 2.
        c_m_negative = np.conjugate(c_m_positive[1:])
        c_m = np.hstack((c_m_negative[::-1], c_m_positive)) / float(self.N)

        if C_window != None:
            # Window the Fourier coefficients.
            # This will damp the highest frequencies

            if self.verbose:
                print('windowing the Fourier coefficients')
            c_m = c_m * c_window(self.m, int(C_window * self.N // 2.))

        A_out = np.zeros((pf.shape[0], self.k_size))
        for i in range(pf.shape[0]):
            # convolve f_c and g_c
            # C_l=np.convolve(c_m*self.g_m[i,:],c_m*self.g_n[i,:])
            C_l = fftconvolve(c_m * g_m[i, :], c_m * g_n[i, :])

            # multiply all l terms together
            C_l = C_l * h_l[i, :] * two_part_l[i]

            # set up to feed ifft an array ordered with l=0,1,...,-1,...,N/2-1
            c_plus = C_l[self.l >= 0]
            c_minus = C_l[self.l < 0]

            C_l = np.hstack((c_plus[:-1], c_minus))
            A_k = ifft(C_l) * C_l.size  # multiply by size to get rid of the normalization in ifft

            A_out[i, :] = np.real(A_k[::2]) * pf[i] * self.k ** (-p[i] - 2)
        # note that you have to take every other element
        # in A_k, due to the extended array created from the
        # discrete convolution

        P_out = irfft(c_m[self.m >= 0]) * self.k ** nu * float(self.N)
        if self.n_pad is not None:
            # get rid of the elements created from padding
            P_out = P_out[self.id_pad]
            A_out = A_out[:, self.id_pad]

        return P_out, A_out

    def J_k_tensor(self, P, X, P_window=None, C_window=None):

        pf, p, nu1, nu2, g_m, g_n, h_l = X

        if self.low_extrap is not None:
            P = self.EK.extrap_P_low(P)

        if self.high_extrap is not None:
            P = self.EK.extrap_P_high(P)

        A_out = np.zeros((pf.size, self.k_size))

        P_fin = np.zeros(self.k_size)

        for i in range(pf.size):

            P_b1 = P * self.k_old ** (-nu1[i])
            P_b2 = P * self.k_old ** (-nu2[i])

            if P_window != None:
                # window the input power spectrum, so that at high and low k
                # the signal smoothly tappers to zero. This make the input
                # more "like" a periodic signal

                if self.verbose:
                    print('windowing biased power spectrum')
                W = p_window(self.k_old, P_window[0], P_window[1])
                P_b1 = P_b1 * W
                P_b2 = P_b2 * W

            if self.n_pad != 0:
                P_b1 = np.pad(P_b1, pad_width=(self.n_pad, self.n_pad), mode='constant', constant_values=0)
                P_b2 = np.pad(P_b2, pad_width=(self.n_pad, self.n_pad), mode='constant', constant_values=0)
            c_m_positive = rfft(P_b1)
            c_n_positive = rfft(P_b2)

            c_m_negative = np.conjugate(c_m_positive[1:])
            c_n_negative = np.conjugate(c_n_positive[1:])
            c_m = np.hstack((c_m_negative[::-1], c_m_positive)) / float(self.N)
            c_n = np.hstack((c_n_negative[::-1], c_n_positive)) / float(self.N)

            if C_window != None:
                # window the Fourier coefficients.
                # This will damping the highest frequencies
                if self.verbose:
                    print('windowing the Fourier coefficients')
                c_m = c_m * c_window(self.m, int(C_window * self.N / 2.))
                c_n = c_n * c_window(self.m, int(C_window * self.N / 2.))

            # convolve f_c and g_c
            C_l = fftconvolve(c_m * g_m[i, :], c_n * g_n[i, :])
            # C_l=convolve(c_m*self.g_m[i,:],c_m*self.g_n[i,:])

            # multiply all l terms together
            # C_l=C_l*self.h_l[i,:]*self.two_part_l[i]
            C_l = C_l * h_l[i, :]
            # set up to feed ifft an array ordered with l=0,1,...,-1,...,N/2-1
            c_plus = C_l[self.l >= 0]
            c_minus = C_l[self.l < 0]

            C_l = np.hstack((c_plus[:-1], c_minus))
            A_k = ifft(C_l) * C_l.size  # multiply by size to get rid of the normalization in ifft

            A_out[i, :] = np.real(A_k[::2]) * pf[i] * self.k ** (p[i])
            # note that you have to take every other element
            # in A_k, due to the extended array created from the
            # discrete convolution
            P_fin += A_out[i, :]
        # P_out=irfft(c_m[self.m>=0])*self.k**self.nu*float(self.N)
        if self.n_pad != 0:
            # get rid of the elements created from padding
            # P_out=P_out[self.id_pad]
            A_out = A_out[:, self.id_pad]
            P_fin = P_fin[self.id_pad]

        return P_fin, A_out


### Example script ###
if __name__ == "__main__":
    # An example script to run FASTPT for (P_22 + P_13) and IA.
    # Makes a plot for P_22 + P_13.
    from time import time

    # Version check
    print('This is FAST-PT version', __version__)

    # load the data file

    d = np.loadtxt('Pk_test.dat')
    # declare k and the power spectrum
    k = d[:, 0]
    P = d[:, 1]

    # set the parameters for the power spectrum window and
    # Fourier coefficient window
    # P_window=np.array([.2,.2])
    C_window = .75
    # document this better in the user manual

    # padding length
    n_pad = int(0.5 * len(k))
    #	to_do=['one_loop_dd','IA_tt']
    to_do = ['one_loop_dd']
    #	to_do=['dd_bias','IA_all']
    # to_do=['all']

    # initialize the FASTPT class
    # including extrapolation to higher and lower k
    t1 = time()
    fastpt = FASTPT(k, to_do=to_do, low_extrap=-5, high_extrap=3, n_pad=n_pad)

    t2 = time()
    # calculate 1loop SPT (and time the operation)
    # P_spt=fastpt.one_loop_dd(P,C_window=C_window)
    P_lpt = fastpt.one_loop_dd_bias_lpt(P, C_window=C_window)

    # for M = 10**14 M_sun/h
    b1L = 1.02817
    b2L = -0.0426292
    b3L = -2.55751
    b1E = 1 + b1L

    # for M = 10**14 M_sun/h
    # b1_lag = 1.1631
    # b2_lag = 0.1162

    # [Ps, Pnb, Pb1L, Pb1L_2, Pb1L_b2L, Pb2L, Pb2L_2, Pb3L, Pb1L_b3L] = [P_lpt[0],P_lpt[1],P_lpt[2],P_lpt[3],P_lpt[4],P_lpt[5],P_lpt[6],P_lpt[7],P_lpt[8]]
    [Ps, Pnb, Pb1L, Pb1L_2, Pb1L_b2L, Pb2L, Pb2L_2] = [P_lpt[0], P_lpt[1], P_lpt[2], P_lpt[3], P_lpt[4], P_lpt[5],
                                                       P_lpt[6]]

    # Pgg_lpt = (b1E**2)*P + Pnb + (b1L)*(Pb1L) + (b1L**2)*Pb1L_2 + (b1L*b2L)*Pb1L_b2L + (b2L)*(Pb2L) + (b2L**2)*Pb2L_2 + (b3L)*(Pb3L) + (b1L*b3L)*Pb1L_b3L
    Pgg_lpt = (b1E ** 2) * P + Pnb + b1L * Pb1L + (b1L ** 2) * Pb1L_2 + (b1L * b2L) * Pb1L_b2L + b2L * Pb2L + (
            b2L ** 2) * Pb2L_2

    # print([pnb,pb1L,pb1L_2,pb2L,Pb1L_b2L])

    t3 = time()
    print('initialization time for', to_do, "%10.3f" % (t2 - t1), 's')
    print('one_loop_dd recurring time', "%10.3f" % (t3 - t2), 's')

    # calculate tidal torque EE and BB P(k)
    # P_IA_tt=fastpt.IA_tt(P,C_window=C_window)
    # P_IA_ta=fastpt.IA_ta(P,C_window=C_window)
    # P_IA_mix=fastpt.IA_mix(P,C_window=C_window)
    # P_RSD=fastpt.RSD_components(P,1.0,C_window=C_window)
    # P_kPol=fastpt.kPol(P,C_window=C_window)
    # P_OV=fastpt.OV(P,C_window=C_window)
    # P_IRres=fastpt.IRres(P,C_window=C_window)
    # make a plot of 1loop SPT results
    import matplotlib.pyplot as plt

    ax = plt.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$P(k)$', size=30)
    ax.set_xlabel(r'$k$', size=30)

    ax.plot(k, P, label='linear')
    # ax.plot(k,P_spt[0], label=r'$P_{22}(k) + P_{13}(k)$' )
    ax.plot(k, Pgg_lpt, label='P_lpt')

    plt.legend(loc=3)
    plt.grid()
    plt.show()
