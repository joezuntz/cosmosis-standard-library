from cosmosis.datablock import option_section, names


def setup(options):
    sections = [
        "shear_calibration_parameters",
        "shear_cl",
        "shear_cl_gg",
        "shear_cl_gi",
        "shear_xi",
        "shear_cl_ii",
        "wl_number_density",
        "wl_photoz_errors",
    ]
    return sections


def execute(block, config):
    sections = config
    for s in sections:
        block._delete_section(s)
    like = block[names.likelihoods, "xipm_like"]
    block._delete_section(names.likelihoods)
    block[names.likelihoods, "xipm_1_like"] = like
    return 0


def cleanup(config):
    pass
