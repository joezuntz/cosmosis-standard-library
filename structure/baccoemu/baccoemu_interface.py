from cosmosis.datablock import option_section, names
from scipy.interpolate import interp2d
import numpy as np
import baccoemu_vendored as baccoemu
print(baccoemu.__file__)


def setup(options):
    mode = options.get_string(option_section, "mode", default="pk_nl_only")

    # do this only once in the setup phase
    emulator = baccoemu.Matter_powerspectrum()

    allowed_modes = ["pk_nl_only", "pk_nl_plus_baryons"]
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
        if mode == "pk_nl_only":
            run_emulator_no_baryons(block, emulator)
        elif mode == "pk_nl_plus_baryons":
            run_emulator_with_baryons(block, emulator)
        else:
            raise RuntimeError("This should not happen")
    except ValueError as error:
        print(error)
        return 1

    return 0


def run_emulator_with_baryons(block, emulator):
    P_lin = block["matter_power_lin", "p_k"]
    k_lin = block["matter_power_lin", "k_h"]
    z_lin = block["matter_power_lin", "z"]

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

    I_nl = interp2d(np.log10(k_bacco), z_lin[zmask], np.log10(F))
    F_interp = 10 ** I_nl(np.log10(k_nl), z_nl)

    # same thing for baryons
    k_bacco, S = emulator.get_baryonic_boost(k=k_lin[kmask], **params)

    I_baryon = interp2d(np.log10(k_bacco), z_lin[zmask], S)
    S_interp = I_baryon(np.log10(k_nl), z_nl)

    # apply the factor
    P_nl = F_interp * S_interp * P_lin

    block._copy_section("matter_power_lin", "matter_power_nl")

    # save everything
    block["matter_power_nl", "p_k"] = P_nl
    block["matter_power_nl", "k_h"] = k_nl
    block["matter_power_nl", "k"] = z_nl

    return 0


def run_emulator_no_baryons(block, emulator):
    P_lin = block["matter_power_lin", "p_k"]
    k_lin = block["matter_power_lin", "k_h"]
    z_lin = block["matter_power_lin", "z"]

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
    I = interp2d(np.log10(k_bacco), z_lin[zmask], np.log10(F))
    F_interp = 10 ** I(np.log10(k_nl), z_nl)

    # apply the factor
    P_nl = F_interp * P_lin

    block._copy_section("matter_power_lin", "matter_power_nl")

    # save everything
    block["matter_power_nl", "p_k"] = P_nl
    block["matter_power_nl", "k_h"] = k_nl
    block["matter_power_nl", "k"] = z_nl

    return 0
