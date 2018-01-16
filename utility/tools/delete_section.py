from __future__ import print_function
from cosmosis.datablock import option_section


def setup(options):
    sections = options.get_string(option_section, "sections")
    sections = sections.split()
    if not sections:
        print("WARNING: No sections specified to delete in delete_section")
    return sections


def execute(block, config):
    sections = config
    for s in sections:
        block._delete_section(s)
    return 0


def cleanup(config):
    pass
