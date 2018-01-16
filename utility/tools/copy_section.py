from builtins import str
from cosmosis.datablock import option_section


def setup(options):
    source = str(options[option_section, "source"])
    dest = str(options[option_section, "dest"])
    return (source, dest)


def execute(block, config):
    source, dest = config
    block._copy_section(source, dest)
    return 0


def cleanup(config):
    pass
