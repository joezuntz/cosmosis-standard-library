from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
    # Get the name of the bias section
    bias_section = options.get_string(option_section, "bias_section", "bias_lens")
	# Get the name of the xlens section
    r_section = options.get_string(option_section,"r_section", "x_lens")

    return bias_section, r_section