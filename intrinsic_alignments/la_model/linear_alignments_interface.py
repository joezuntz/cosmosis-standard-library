from cosmosis.datablock import names, option_section
from linear_alignments import kirk_rassat_host_bridle_power
from linear_alignments import bridle_king
from linear_alignments import bridle_king_corrected
from linear_alignments import linear
import numpy as np


def setup(options):
    method = options[option_section, "method"].lower()
    grid_mode = options.get_bool(option_section, "grid_mode", default=False)
    gal_intrinsic_power = options.get_bool(
        option_section, "do_galaxy_intrinsic", False)
    name = options.get_string(option_section, "name", default="").lower()
    if name:
        suffix = "_" + name
    else:
        suffix = ""

    if method not in ["krhb", "bk", "bk_corrected", "linear"]:
        raise ValueError('The method in the linear alignment module must'
                         'be either "KRHB" (for Kirk, Rassat, Host, Bridle) or BK for '
                         'Bridle & King or "BK_corrected" for the corrected version of that')

    return method, gal_intrinsic_power, grid_mode, suffix


def execute(block, config):
    method, gal_intrinsic_power, grid_mode, suffix = config

    # load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
    lin = names.matter_power_lin
    nl = names.matter_power_nl
    ia = names.intrinsic_alignment_parameters + suffix
    ia_ii = names.intrinsic_power + suffix
    ia_gi = names.galaxy_intrinsic_power + suffix
    ia_mi = names.matter_intrinsic_power + suffix
    gm = "matter_galaxy_power" + suffix
    cosmo = names.cosmological_parameters

    z_lin, k_lin, p_lin = block.get_grid(lin, "z", "k_h", "p_k")
    z_nl, k_nl, p_nl = block.get_grid(nl, "z", "k_h", "p_k")

    omega_m = block[cosmo, "omega_m"]
    A = block[ia, "A"]

    # run computation and write to datablock
    if method == 'krhb':
        P_II, P_GI, b_I, r_I, k_I, z_I = kirk_rassat_host_bridle_power(
            z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)
    elif method == 'bk':
        P_II, P_GI, b_I, r_I, k_I, z_I = bridle_king(
            z_nl, k_nl, p_nl, A, omega_m)
    elif method == 'bk_corrected':
        P_II, P_GI, b_I, r_I, k_I, z_I = bridle_king_corrected(
            z_nl, k_nl, p_nl, A, omega_m)
    elif method == "linear":
        P_II, P_GI, b_I, r_I, k_I, z_I = linear(
            z_lin, k_lin, p_lin, A, omega_m)

    if grid_mode:
        block.put_grid(ia, "z", z_I, "k_h", k_I,  "b_I" + suffix, b_I)
        block.put_grid(ia, "z", z_I, "k_h", k_I, "r_I" + suffix, r_I)
    else:
        block.put_grid(ia_mi, "z", z_I, "k_h", k_I,  "p_k", P_GI)
        block.put_grid(ia_ii, "z", z_I, "k_h", k_I, "p_k", P_II)

    # This is a bit of hack...scale GI power spectrum (which is really matter-intrinsic
    # power spectrum) by P_gal_matter/P_delta_delta
    if gal_intrinsic_power:
        z, k, p_gm = block.get_grid(gm, "z", "k_h", "p_k")
        P_gI = P_GI * p_gm / p_nl
        block.put_grid(ia_gi, "z", z, "k_h", k, "p_k", P_gI)
    return 0


def cleanup(config):
    pass
