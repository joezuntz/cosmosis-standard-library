from cosmosis.datablock import names

halo_model_params = names.halo_model_parameters


def setup(options):
    return {}

def execute(block, config):
    A = block[halo_model_params, "A"]
    eta0 = 0.98 - 0.12 * A
    block[halo_model_params, "eta"] = eta0
    return 0
