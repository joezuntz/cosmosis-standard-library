import cosmosis_py.section_names
from cosmosis_py.block import option_section
import cfhtlens_like
from cfhtlens_like import n_z_bin

def setup(options):
    sec = option_section
    covmat_file = options.get(sec, 'covariance_file', default=cfhtlens_like.DEFAULT_COVMAT)
    data_file = options.get(sec, 'data_file', default=cfhtlens_like.DEFAULT_DATA)
    xiplus_only = options.get(sec, 'xplus_only', default=False)
    cut_low_theta = options.get(sec, 'cut_low_theta', default=True)

    #create likelihood calculator
    #loads named files and prepares itself
    calculator = cfhtlens_like.CFHTLensLikelihood(covmat_file, data_file,plus_only,cut_low_theta)

    #pass back config to the 
    return calculator


def execute(block, config):
    calculator = config

    #Get theta for the sample values
    section=cosmosis_py.section_names.shear_xi
    theta = block[section, "theta"]

    #Get the xi(theta) for these samples, for each pair of bins.
    #The likelihood calculator wants a big dictionary
    xi_data = {}
    for i in xrange(1, n_z_bin+1):
        for j in xrange(i, n_z_bin+1):
            name = 'bin_%d_%d' % (j,i)
            xi = block[section, name]
            xi_data[(i,j)] = (theta, xi)

    #Calculate the likelihood
    like = calculator(xi_data)

    #save the result
    section=pydesglue.section_names.likelihoods
    block[section, "cfhtls_like"] = like

    return 0


