from cosmosis.datablock import option_section, names
from scipy.interpolate import RectBivariateSpline
import numpy as np
import baccoemu_vendored as baccoemu
import traceback


def setup(options):
    mode = options.get_string(option_section, "mode", default="nonlinear")

    # do this only once in the setup phase
    emulator = baccoemu.Matter_powerspectrum()

    allowed_modes = ["nonlinear", "baryons", "nonlinear+baryons"]
    if mode not in allowed_modes:
        raise ValueError(f"unrecognised value of 'mode' parameter in baccoemu: {mode}")

    return mode, emulator


def check_params(params):
    ranges = {
        "omega_cold": [0.15, 0.6],
        "omega_baryon": [0.03, 0.07],
        "ns": [0.92, 1.01],
        "hubble": [0.6, 0.8],
        "neutrino_mass": [0.0, 0.4],
    }

    for key, (vmin, vmax) in ranges.items():
        if not (params[key] > vmin) and (params[key] < vmax):
            raise ValueError(f"BaccoEmu: {key} must be in range [{vmin},{vmax}]")


def execute(block, config):
    mode, emulator = config

    try:
        if mode == "nonlinear":
            emulate_nonlinear_power(block, emulator)
        elif mode == "nonlinear+baryons":
            emulate_boosted_nonlinear_power(block, emulator)
        elif mode == "baryons":
            emulate_baryonic_boost(block, emulator)
        else:
            raise RuntimeError("This should not happen")
    except ValueError as error:
        # print traceback from exception
        print(traceback.format_exc())
        return 1

    return 0


def emulate_nonlinear_power(block, emulator):
    z_lin, k_lin, P_lin = block.get_grid("matter_power_lin", "z", "k_h", "p_k")

    # bacco specific stuff
    # This is required to avoid a bounds error
    zmask = z_lin < 1.5
    a_lin = 1 / (1 + z_lin[zmask])

    kmask = (k_lin < 4.69) & (k_lin > 0.0001)
    cosmo = names.cosmological_parameters

    omega_cold = block[cosmo, "omega_c"] + block[cosmo, "omega_b"]
    Omb = block[cosmo, "omega_b"]

    params = {
        "omega_cold": omega_cold,
        "A_s": block[cosmo, "a_s"],
        "omega_baryon": Omb,
        "ns": block[cosmo, "n_s"],
        "hubble": block[cosmo, "h0"],
        "neutrino_mass": block[cosmo, "mnu"],
        "w0": block.get_double(cosmo, "w", -1.0),
        "wa": 0.0,
        "expfactor": a_lin,
    }

    # check we're within the allowed bounds for the emulator
    check_params(params)

    # evaluate the nonlinear growth factor as a function of k and z
    k_bacco, F = emulator.get_nonlinear_boost(k=k_lin[kmask], cold=False, **params)
    k_nl = k_lin
    z_nl = z_lin

    # interpolate it to the same sampling as the linear matter power spectrum
    I = RectBivariateSpline(np.log10(k_bacco), z_lin[zmask], np.log10(F).T)
    F_interp = 10 ** I(np.log10(k_nl), z_nl).T

    # apply the factor
    P_nl = F_interp * P_lin

    # save everything
    block.put_grid("matter_power_nl", "z", z_nl, "k_h", k_nl, "p_k", P_nl)

    return 0


def emulate_boosted_nonlinear_power(block, emulator):
    z_lin, k_lin, P_lin = block.get_grid("matter_power_lin", "z", "k_h", "p_k")

    # This is required to avoid a bounds error
    zmask = z_lin < 1.5
    a_lin = 1 / (1 + z_lin[zmask])

    kmask = (k_lin < 4.69) & (k_lin > 0.0001)

    cosmo = names.cosmological_parameters

    # In BACCO omega_cold refers to baryons + CDM
    omega_cold = block[cosmo, "omega_c"] + block[cosmo, "omega_b"]

    params = {
        "omega_cold": omega_cold,
        "A_s": block[cosmo, "a_s"],
        "omega_baryon": block[cosmo, "omega_b"],
        "ns": block[cosmo, "n_s"],
        "hubble": block[cosmo, "h0"],
        "neutrino_mass": block[cosmo, "mnu"],
        "w0": block.get_double(cosmo, "w", -1.0),
        "wa": block.get_double(cosmo, "wa", 0.0),
        "expfactor": a_lin,
        "M_c": block["baryon_parameters", "m_c"],
        "eta": block["baryon_parameters", "eta"],
        "beta": block["baryon_parameters", "beta"],
        "M1_z0_cen": block["baryon_parameters", "m1_z0_cen"],
        "theta_out": block["baryon_parameters", "theta_out"],
        "theta_inn": block["baryon_parameters", "theta_inn"],
        "M_inn": block["baryon_parameters", "m_inn"],
    }

    # check we're within the allowed bounds for the emulator
    check_params(params)

    k_nl = k_lin
    z_nl = z_lin

    # evaluate the nonlinear growth factor as a function of k and z
    k_bacco, F = emulator.get_nonlinear_boost(k=k_lin[kmask], cold=False, **params)
    I_nl = RectBivariateSpline(np.log10(k_bacco), z_lin[zmask], np.log10(F).T)
    F_interp = 10 ** I_nl(np.log10(k_nl), z_nl).T

    # same thing for baryons
    k_bacco, S = emulator.get_baryonic_boost(k=k_lin[kmask], **params)
    I_baryon = RectBivariateSpline(np.log10(k_bacco), z_lin[zmask], S.T)
    S_interp = I_baryon(np.log10(k_nl), z_nl).T

    # apply the factor
    P_nl = F_interp * S_interp * P_lin

    # save everything
    block.put_grid("matter_power_nl", "z", z_nl, "k_h", k_nl, "p_k", P_nl)

    return 0


def emulate_baryonic_boost(block, emulator):
    z_nl, k_nl, P_nl = block.get_grid("matter_power_nl", "z", "k_h", "p_k")

    # This is required to avoid a bounds error
    zmask = z_nl < 1.5
    a_lin = 1 / (1 + z_nl[zmask])

    kmask = (k_nl < 4.69) & (k_nl > 0.0001)

    cosmo = names.cosmological_parameters

    # In BACCO omega_cold refers to baryons + CDM
    omega_cold = block[cosmo, "omega_c"] + block[cosmo, "omega_b"]

    params = {
        "omega_cold": omega_cold,
        "A_s": block[cosmo, "a_s"],
        "omega_baryon": block[cosmo, "omega_b"],
        "ns": block[cosmo, "n_s"],
        "hubble": block[cosmo, "h0"],
        "neutrino_mass": block[cosmo, "mnu"],
        "w0": block.get_double(cosmo, "w", -1.0),
        "wa": block.get_double(cosmo, "wa", 0.0),
        "expfactor": a_lin,
        "M_c": block["baryon_parameters", "m_c"],
        "eta": block["baryon_parameters", "eta"],
        "beta": block["baryon_parameters", "beta"],
        "M1_z0_cen": block["baryon_parameters", "m1_z0_cen"],
        "theta_out": block["baryon_parameters", "theta_out"],
        "theta_inn": block["baryon_parameters", "theta_inn"],
        "M_inn": block["baryon_parameters", "m_inn"],
    }

    # check we're within the allowed bounds for the emulator
    check_params(params)

    # just get the boost factor from bacco
    k_bacco, S = emulator.get_baryonic_boost(k=k_nl[kmask], **params)
    I_baryon = RectBivariateSpline(np.log10(k_bacco), z_nl[zmask], S.T)
    S_interp = I_baryon(np.log10(k_nl), z_nl).T

    # apply the factor
    P_nl *= S_interp

    # save everything
    block.put_grid("matter_power_nl", "z", z_nl, "k_h", k_nl, "p_k", P_nl)

    return 0
