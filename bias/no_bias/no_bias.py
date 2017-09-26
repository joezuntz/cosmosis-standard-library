from __future__ import print_function
from cosmosis.datablock import names, option_section


def setup(options):
    use_lin_power = options.get_bool(option_section, "use_lin_power", False)
    return use_lin_power


def execute(block, config):
    use_lin_power = config
    if use_lin_power:
        pk_section = names.matter_power_lin
    else:
        pk_section = names.matter_power_nl
    block._copy_section(pk_section, names.galaxy_power)
    # if using use_lin_power, we're assuming galaxies follow linear fluctuations
    # then the matter-galaxy cross power should be sqrt(P_lin*P_nl)
    if use_lin_power:
        print('WARNING: galaxy power and matter-galaxy power no consistent')
    block._copy_section(names.matter_power_nl, names.matter_galaxy_power)
    if block.has_section("matter_intrinsic_power"):
        block._copy_section("matter_intrinsic_power", "galaxy_intrinsic_power")
    return 0
