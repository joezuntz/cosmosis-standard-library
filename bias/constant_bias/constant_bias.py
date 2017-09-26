from cosmosis.datablock import names, option_section


def setup(options):
    return []


def execute(block, config):
    b = block['galaxy_bias', 'b']

    # Copy the non-linear matter power to both galaxy-power and
    # matter-galaxy cross power (the latter is used in lensing-position spectra)
    block._copy_section(names.matter_power_nl, names.galaxy_power)
    block._copy_section(names.matter_power_nl, names.matter_galaxy_power)

    # Now apply constant biases to the values we have just copied.
    # More realistic bias models can use b(k,z) for example.
    # The cross power of course picks up only one factor of bias.
    block[names.galaxy_power, "P_K"] *= b**2
    block[names.matter_galaxy_power, "P_K"] *= b

    # We may have a matter intrinsic power aleady worked out.
    # Copy that if so.
    if block.has_section("matter_intrinsic_power"):
        block._copy_section("matter_intrinsic_power", "galaxy_intrinsic_power")

        # and apply the bias here too, again only one factor.
        block[names.galaxy_intrinsic_power, "P_K"] *= b

    return 0
