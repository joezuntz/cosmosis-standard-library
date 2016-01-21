from cosmosis.datablock import names, option_section

def setup(options):
    return []

def execute(block, config):
    block._copy_section(names.matter_power_nl, names.galaxy_power)
    block._copy_section(names.matter_power_nl, names.matter_galaxy_power)
    return 0