from cosmosis.datablock import option_section
import re


def setup(options):
    source = str(options[option_section, "source"])
    dest = str(options[option_section, "dest"])

    # Split on either comma or space
    source = re.split("[ ,]+", source)
    dest = re.split("[ ,]+", dest)
    if len(source) != len(dest):
        raise ValueError("In the rename module the source and dest parameters must "
                         "have the same number of sections in")
    return (source, dest)


def execute(block, config):
    source, dest = config
    for (s, d) in zip(source, dest):
        block._copy_section(s, d)
        block._delete_section(s)
    return 0


def cleanup(config):
    pass
