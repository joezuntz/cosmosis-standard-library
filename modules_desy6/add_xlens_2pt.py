from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
	# Get the name of the xlens section
    x_section = options.get_string(option_section,"xlens", "xlens")

    return x_section


def execute(block, config):

    x_section = config

    xlens_val = block[x_section, 'xlens_all']


    for i in range(1,6):
        for j in range(1,4):
            if j >= i:
                ggl_bin_label = 'bin_'+str(i)+'_'+str(j)
                ggl = block['galaxy_shear_xi',ggl_bin_label]
                ggl_xlens = xlens_val*ggl
                block['galaxy_shear_xi', ggl_bin_label] = ggl_xlens

    return 0
