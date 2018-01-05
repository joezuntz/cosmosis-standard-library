from builtins import zip
from builtins import str
from cosmosis.datablock import option_section


def setup(options):
    source = str(options[option_section, "source"])
    source = source.split(',')
    dest = str(options[option_section, "dest"])
    dest = dest.split(',')
    return (source, dest)


def execute(block, config):
    source, dest = config
    for (s, d) in zip(source, dest):
        block._copy_section(s, d)
        block._delete_section(s)
    return 0


def cleanup(config):
    pass
