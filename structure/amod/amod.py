from cosmosis.datablock import option_section, names


def setup(options):
    return None


def execute(block):
    # inputs
    pk_lin = block["matter_power_lin", "p_k"]
    pk_nonlin = block["matter_power_nl", "p_k"]
    A_mod = block["amod_parameter", "a_mod"]

    # for now assume pk_lin and pk_nonlin are defined on the same grids
    assert (block["matter_power_lin", "z"] == block["matter_power_nl", "z"]).all()
    assert (block["matter_power_lin", "k_h"] == block["matter_power_nl", "k_h"]).all()

    # calculate new non-linear power
    pk_m = pk_lin + A_mod * (pk_nonlin - pk_lin)
    block["matter_power_nl", "p_k"] = pk_m

    return 0


def cleanup(config):
    pass
