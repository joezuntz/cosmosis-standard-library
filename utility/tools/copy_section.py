from builtins import str
from cosmosis.datablock import option_section


def setup(options):
    sources = str(options[option_section, "source"]).split()
    dests = str(options[option_section, "dest"]).split()
    for source, dest in zip(sources, dests):
        print("Will copy section {} -> {}".format(source, dest))
    return (sources, dests)


def execute(block, config):
    sources, dests = config
    for source, dest in zip(sources, dests):
        block._copy_section(source, dest)
    return 0


def cleanup(config):
    pass
