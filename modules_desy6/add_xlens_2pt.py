from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
	# Get the name of the xlens section
    r_section = options.get_string(option_section,"r_section", "x_lens")

    return bias_section, r_section


def execute(block, config):

    r_section = config

    rval = block[r_section, 'rmean_bin']


    for i in range(1,6):
        for j in range(1,6):
            if j >= i:
                ggl_bin_label = 'bin_'+str(i)+'_'+str(j)
                ggl = block['galaxy_shear_xi',ggl_bin_label]
                ggl_xlens = rval*ggl
                block['galaxy_shear_xi', ggl_bin_label] = ggl_xlens

    return 0
