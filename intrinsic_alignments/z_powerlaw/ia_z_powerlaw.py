from cosmosis.datablock import names, option_section
import numpy as np


def setup(options):
    name = options.get_string(option_section, "name", default="").lower()
    do_galaxy_intrinsic = options.get_bool(
        option_section, "do_galaxy_intrinsic", default=False)
    if name:
        suffix = "_" + name
    else:
        suffix = ""
    return {"suffix": suffix, "do_galaxy_intrinsic": do_galaxy_intrinsic}


def execute(block, config):
    # Get the names of the sections to save to
    suffix = config['suffix']
    ia_section = names.intrinsic_alignment_parameters + suffix
    ia_ii = names.intrinsic_power + suffix
    ia_mi = names.matter_intrinsic_power + suffix

    # read in power spectra for P_II and P_MI
    z, k, p_ii = block.get_grid(ia_ii, "z", "k_h", "p_k")
    z, k, p_mi = block.get_grid(ia_mi, "z", "k_h", "p_k")

    # read alpha from ia_section values section
    alpha = block[ia_section, 'alpha']
    z0 = block.get_double(ia_section, 'z0', default=0.0)
    _, z_grid = np.meshgrid(k, z)

    # Construct and apply redshift scaling
    z_scaling = ((1 + z_grid) / (1 + z0))**alpha
    p_ii *= z_scaling**2
    p_mi *= z_scaling

    # Save grid back to the block
    block.replace_grid(ia_ii, "z", z, "k_h", k, "p_k", p_ii)
    block.replace_grid(ia_mi, "z", z, "k_h", k, "p_k", p_mi)

    if config['do_galaxy_intrinsic']:
        ia_gi = names.galaxy_intrinsic_power + suffix
        z, k, p_gi = block.get_grid(ia_gi, "z", "k_h", "p_k")
        p_gi *= z_scaling

        block.replace_grid(ia_gi, "z", z, "k_h", k, "p_k", p_gi)

    return 0


def cleanup(config):
    pass
