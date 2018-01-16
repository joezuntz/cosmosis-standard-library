from __future__ import print_function
from builtins import range
from cosmosis.datablock import option_section, names
import clerkin
import numpy as np

models = ['gtd', 'q', 'q-gtd']  # ...]
modes = ["power", "bias", "both"]


def setup(options):
    model = options.get_string(option_section, "model", default="gtd").lower()
    bias_only = options.get_bool(option_section, "bias_only", default=False)
    suffix = options.get_string(option_section, "suffix", default="")

    if model not in models:
        raise ValueError(
            "The Clerkin module requires a model chosen from: %r" % models)

    if suffix != '':
        suffix = '_' + suffix

    return model, bias_only, suffix


def execute_power(block, model, suffix):

    # Get matter power
    k, z, P = block.get_grid(names.matter_power_nl, "k_h", "z", "p_k")

    if 'gtd' in model:
        # Get growth
        try:
            z_growth = block[names.growth_parameters, 'z']
            growth = block[names.growth_parameters, 'd_z']
        except:
            raise ValueError("The Clerkin module needs you to calculate the growth factor first." +
                             "  Add the growth module before it in your module list")

        # Get bias parameters
        alpha = block['galaxy_bias', 'alpha']
        b0 = block['galaxy_bias', 'b0']
        c = block['galaxy_bias', 'c']

    if 'q' in model:
        # Need amplitude and quadratic part for the scale dependence in this model
        Q = block['galaxy_bias', 'Q']
        A = block['galaxy_bias', 'A']

    # Get matter power with bias
    if model == 'gtd':
        bias = clerkin.gtd_model(k, z, z_growth, growth, alpha, b0, c)
    elif model == 'q':
        bias = clerkin.q_model(k, z, Q, A)
    elif model == 'q-gtd':
        bias = clerkin.gtd_q_model(k, z, z_growth, growth, alpha, b0, c, Q, A)
    else:
        raise ValueError("Unknown model in Clerkin")

    if (bias < 0).any():
        print("Negative galaxy power from bias model (bad model; returning error).")
        return 1

    P_galaxy = bias**2 * P
    P_matter_galaxy = bias * P

    block.put_grid(names.bias_field + suffix, "k_h", k, "z", z, "b", bias)
    block.put_grid("galaxy_power" + suffix, "k_h", k, "z", z, "P_k", P_galaxy)
    block.put_grid("matter_galaxy_power" + suffix, "k_h",
                   k, "z", z, "P_k", P_matter_galaxy)

    return 0


def execute_bias_only(block, model, suffix):
    if model == "gtd":
        # Get growth
        z1 = block[names.growth_parameters, 'z']
        growth_z = block[names.growth_parameters, 'd_z']
        # Get bias parameters
        alpha = block['bias_parameters', 'alpha']
        b0 = block['bias_parameters', 'b0']
        c = block['bias_parameters', 'c']
        b = clerkin.gtd_bias(z1, growth_z, alpha, b0, c)

        nk = len(block['matter_power_nl', 'k_h'])
        bias = np.array([np.array(b) for i in range(nk)])
        #bias = bias.swapaxes(0,1)

        k_h = block['matter_power_nl', 'k_h']
        block.put_grid(names.bias_field, "k_h", k_h,
                       "z", z1, "b" + suffix,  bias)
        block.put_grid(names.bias_field, "k_h", k_h, "z",
                       z1, "r" + suffix, np.ones_like(bias))

        #import pdb ; pdb.set_trace()

    elif model == "q":
        k = np.logspace(-6, 2, 500)
        Q = block['bias_parameters', 'Q']
        A = block['bias_parameters', 'A']
        b = clerkin.q_bias(k, Q, A)
        block[names.bias_field, "k"] = k
        block[names.bias_field, "b" + suffix] = b

    elif model == 'q-gtd':
        z1 = block[names.growth_parameters, 'z']
        growth_z = block[names.growth_parameters, 'd_z']
        # Get bias parameters
        alpha = block['bias_parameters', 'alpha']
        b0 = block['bias_parameters', 'b0']
        c = block['bias_parameters', 'c']
        b = clerkin.gtd_bias(z1, growth_z, alpha, b0, c)
        k = np.logspace(-6, 2, 500)
        Q = block['bias_parameters', 'Q']
        A = block['bias_parameters', 'A']
        b_z = clerkin.gtd_bias(z1, growth_z, alpha, b0, c)
        b_k = clerkin.q_bias(k, Q, A)
        nk = len(k)
        nz = len(z1)
        b = np.zeros((nk, nz))
        for i in range(nk):
            for j in range(nz):
                b[i, j] = b_k[i] * b_z[j]

        block.put_grid(names.bias_field, "k_h", k, "z", z1, "b" + suffix, b)
    return 0


def execute_both(block, model, suffix):
    status = execute_bias(block, model, suffix)
    if status:
        return status
    return execute_power(block, model, suffix)


def execute(block, config):
    model, bias_only, suffix = config
    if bias_only:
        return execute_bias_only(block, model, suffix)
    else:
        return execute_power(block, model, suffix)
