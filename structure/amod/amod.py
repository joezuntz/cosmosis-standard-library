from cosmosis.datablock import option_section, names


def setup(options):
    return None


def execute(block):
    # inputs
    pk_lin = block["matter_power_lin", "p_k"]
    pk_nonlin = block["matter_power_nl", "p_k"]
    A_mod = block["a_mod_parameter", "a_mod"]

    # calculate new non-linear power
    pk_m = pk_lin + A_mod * (pk_nonlin - pk_lin)
    block["matter_power_nl", "p_k"] = pk_m

    return 0


def cleanup(config):
    pass
