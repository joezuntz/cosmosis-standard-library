"""
An interface to the determine the SP(k) model of baryon feedback 
suppression to the non-linear matter power spectrum.

"""

from cosmosis.datablock import option_section, names as section_names
from cosmosis.utils import datablock_to_astropy  #new for cosmosis v3.13

import numpy as np
import astropy.cosmology
import warnings
from scipy.interpolate import LinearNDInterpolator

# pyspk tells the whole of python to always print every warning over and over again.
# Undo it like this so the warning is only printed out once
warnings_filters = warnings.filters[:]
import pyspk.model as spk
warnings.filters = warnings_filters


def my_ceil(a, precision=0):
    return np.true_divide(np.ceil(a * 10**precision), 10**precision)


def my_floor(a, precision=0):
    return np.true_divide(np.floor(a * 10**precision), 10**precision)


round_up = np.vectorize(my_ceil)
round_down = np.vectorize(my_floor)


def setup(options):

    verbose = options.get_bool(option_section, "verbose", default=False)
    SO = options.get_int(option_section, "SO", default=500)
    fb_table = options.get_string(option_section, "fb_table", default="")

    # If a baryon fraction - halo mass table is provided, load it and create an interpolator
    if fb_table:
        tab = np.loadtxt(fb_table, skiprows=1, delimiter=",")
        inter = LinearNDInterpolator(tab[:, [0, 1]], tab[:, 2], rescale=True)
    else:
        inter = None
        fb_table = None

    extrapolate = options.get_bool(option_section, "extrapolate", default=False)

    # These are parameters for spk
    config = {}
    config["verbose"] = verbose
    config["SO"] = SO
    config["fb_table"] = fb_table
    config["extrapolate"] = extrapolate
    config["interpolator"] = inter

    # Check that the chosen value for the spherical overdensity is either 200 or 500
    if SO not in [200, 500]:
        raise ValueError(
            f"[SPK]: The spherical overdensity SO must be either 200 or 500. You provided: {SO}"
        )

    return config


def check_parameter_choice(fb_table, spk_params):
    """
    Look at the parameters provided for SPK and check that
    they are consistent with exactly one method.
    """
    #SP(k) allows for four different methods to fix the fb-Mhalo relation
    # 1. A power law with 3 parameters, fb_a, fb_pow, fb_pivot and an optional 4th parameter, m_pivot (default 1 M_odot)
    # 2. A redshift dependent power law with 3 parameters: alpha, beta, gamma (here m_pivot is fixed to 10^14 M_odot)
    # 3. A redshift dependent double power law with 5 parameters: epsilon, alpha, beta, gamma, m_pivot
    # 4. A non-parametric table of fb values as a function of redshift and halo mass

    # If fb_table is set, we use method 4 irrespective of the other parameters
    # Warn the user that their values are not being used in this instance
    if fb_table:
        for key, value in spk_params.items():
            if value:
                spk_params[key] = None # Reset the value to None
                warnings.warn(
                    f"[SPK]: As you have set fb_table in the param file, "+key+" will be ignored."
                )
    else:   
        # If m_pivot is set, we use method 1 or 3 depending on the other parameters       
        if spk_params["m_pivot"]:
            modes = [
            (
                "Power law using m_pivot",
                ["fb_a", "fb_pow", "fb_pivot"],
                spk_params["fb_a"],
                spk_params["fb_pow"],
                spk_params["fb_pivot"],
            ),
            (
                "Redshift dependent double power law",
                ["epsilon", "alpha", "beta", "gamma", "m_pivot"],
                spk_params["epsilon"],
                spk_params["alpha"],
                spk_params["beta"],
                spk_params["gamma"],
                spk_params["m_pivot"],
            ),
            ]
        # Otherwise it's methods 1 or 2
        else:
            modes = [
            (
                "Power law using default m_pivot",
                ["fb_a", "fb_pow", "fb_pivot"],
                spk_params["fb_a"],
                spk_params["fb_pow"],
                spk_params["fb_pivot"],
            ),
            (
                "Redshift dependent power law using fixed m_pivot",
                ["alpha", "beta", "gamma"],
                spk_params["alpha"],
                spk_params["beta"],
                spk_params["gamma"],
            ),
            ]

        provided_modes = sum(any(param is not None for param in mode[2:]) for mode in modes)

        # Complain if more than one mode is consistent with the parameters specfified.
        # This will happen if parameters that should be left unset are supplied.
        if provided_modes != 1:
            provided_params = [
                param
                for mode in modes
                for param in mode[1]
                if mode[2:].count(None) != len(mode[2:])
            ]
            raise ValueError(
                f"[SPK]: Only one method to specify the baryon fraction should be provided. You provided: {provided_params}. \
                Please either provide f_a, f_b, f_pivot (with or without m_pivot) or epsilon, alpha, beta, gamma and m_pivot."
            )

        # Now complain if not all parameters are set for a chosen mode.
        for mode in modes:
            if any(param is not None for param in mode[2:]):
                missing_params = [
                    param_name
                    for param_name, param_value in zip(mode[1], mode[2:])
                    if param_value is None
                ]
                if len(missing_params) != 0:
                    raise ValueError(
                        f"[SPK]: The following parameter(s) is(are) missing in method '{mode[0]}': {', '.join(missing_params)}."
                    )

def execute(block, config):
    verbose = config["verbose"]
    SO = config["SO"]
    fb_table = config["fb_table"]
    extrapolate = config["extrapolate"]
    inter = config["interpolator"]

    # Pull the cosmological model from the block in astropy format
    cosmo = datablock_to_astropy(block)

    M_halo = None
    fb = None

    # Load the current non-linear matter power spectrum
    section = section_names.matter_power_nl
    k, z_array, P = block.get_grid(section, "k_h", "z", "P_k")

    params = [
        "fb_a",
        "fb_pow",
        "fb_pivot",
        "epsilon",
        "alpha",
        "beta",
        "gamma",
        "m_pivot",
    ]

    # Get any of the named parameters from above that are in the block
    spk_params = {}
    for param_i in params:
        if block.has_value("spk", param_i):
            spk_params[param_i] = block.get("spk", param_i)
        else:
            spk_params[param_i] = None

    # Check that the specified parameters match exactly one
    # sensible mode.
    check_parameter_choice(fb_table, spk_params)

    sup_array = np.empty_like(P)

    for i in range(len(z_array)):
        z = z_array[i]

        if inter is not None:
            min_mass = round_down(np.log10(spk.optimal_mass(SO, z, np.max(k))), 1)
            max_mass = round_up(np.log10(spk.optimal_mass(SO, z, np.min(k))), 1)
            M_halo = np.logspace(min_mass, max_mass, 100)
            fb = inter(z, M_halo)
            if np.isnan(np.sum(fb)):
                raise ValueError(
                    f"[SPK]: Requested values (z and/or M_halo) outside of the convex hull of the input points in fb_table. Hint: Check your fb_table and the requested redshifts."
                )

        # Get the suppression factor for this redshift
        k, sup = spk.sup_model(
            SO=SO,
            z=z,
            fb_a=spk_params["fb_a"],
            fb_pow=spk_params["fb_pow"],
            fb_pivot=spk_params["fb_pivot"],
            M_halo=M_halo,
            fb=fb,
            extrapolate=extrapolate,
            epsilon=spk_params["epsilon"],
            alpha=spk_params["alpha"],
            beta=spk_params["beta"],
            gamma=spk_params["gamma"],
            m_pivot=spk_params["m_pivot"],
            cosmo=cosmo,
            k_array=k,
            verbose=verbose,
        )
        sup_array[:, i] = sup

    # Scale the matter power spectrum by the suppression
    P_mod = P * sup_array

    # JZ I've also replaced this to avoid assuming the ordering
    # of the P grid in memory.
    section = section_names.matter_power_nl
    block.replace_grid(section, "k_h", k, "z", z_array, "P_k", P_mod)

    # Check for NaNs in the modified power spectrum which occur when 
    # baryon fraction values outside fitting limits are requested
    if np.isnan(P_mod).any(): return 1
    else: return 0
