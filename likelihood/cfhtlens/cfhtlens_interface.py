import cosmosis_py.section_names
from cosmosis_py.block import option_section
import cfhtlens_like
from cfhtlens_like import n_z_bin
import numpy as np

def setup(options):
    sec = option_section
    covmat_file = options.get_string(sec, 'covariance_file', default=cfhtlens_like.DEFAULT_COVMAT)
    data_file = options.get_string(sec, 'data_file', default=cfhtlens_like.DEFAULT_DATA)
    xiplus_only = options.get_bool(sec, 'xplus_only', default=False)
    cut_low_theta = options.get_bool(sec, 'cut_low_theta', default=True)

    #create likelihood calculator
    #loads named files and prepares itself
    calculator = cfhtlens_like.CFHTLensLikelihood(covmat_file, data_file,xiplus_only,cut_low_theta)

    #pass back config to the 
    return calculator


def execute(block, config):
    calculator = config

    #Get theta for the sample values.
    #We need theta twice because our CFHTLens
    #code wants xminus and xplus
    section=cosmosis_py.section_names.shear_xi
    theta = block[section, "theta"]
    theta = np.concatenate((theta, theta))

    #Get the xi(theta) for these samples, for each pair of bins.
    #The likelihood calculator wants a big dictionary
    xi_data = {}
    for i in xrange(1, n_z_bin+1):
        for j in xrange(i, n_z_bin+1):
            name = 'xiplus_%d_%d' % (j,i)
            xiplus = block[section, name]
            name = 'ximinus_%d_%d' % (j,i)
            ximinus = block[section, name]
            xi = np.concatenate((xiplus, ximinus))
            xi_data[(i,j)] = (theta, xi)

    #Calculate the likelihood

    like = calculator(xi_data)

    #save the result
    section=cosmosis_py.section_names.likelihoods
    block[section, "cfhtlens_like"] = like

    return 0



def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness.  The joy of python.
    return 0
