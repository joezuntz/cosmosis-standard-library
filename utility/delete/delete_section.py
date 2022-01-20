from cosmosis.datablock import option_section
import re

def setup(options):
    sections = options.get_string(option_section, "sections")
    # Split on either comma or space
    sections = re.split("[ ,]+", sections)

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
