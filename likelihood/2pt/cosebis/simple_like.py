import numpy as np
from cosmosis.datablock import names
from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors

import twopoint as twopoint

## This simple likelihood calculation is set up for COSEBIs.
## In principle it could be used for any statistic
## provided you don't want to make any cuts to the data, theory and covariance.
## If you want to do something more complex use 2pt_like.py

##--------------------------------------------------
## Function to build the theory Spectrum Measurement 

class SpectrumBuilder():
  
    def __init__(self):
        self.tIL1    = []
        self.tIL2    = []
        self.aIL     = []
        self.angList = []
        self.valList = []
        return
    
    def addTomo(self, tomoInd1, tomoInd2, angle, value):
        N_ang  = len(angle)
        angInd = np.arange(N_ang, dtype=int)
        self.tIL1.append([tomoInd1+1] * N_ang)
        self.tIL2.append([tomoInd2+1] * N_ang)
        self.aIL.append(angInd + 1)
        self.angList.append(angle)
        self.valList.append(value)
        return

##--------------------------------------------------

def setup(options):
    config = {}
    
    ## Read inputs from datablock and ensure data file exists
    try:
        config['data_file'] = options.get_string(option_section, 'data_file')
    except errors.BlockNameNotFound:
        raise NameError('data_file cannot be empty')

    ## Read data choice names for data inputs
    config['stat_name'], config['angle_name']= options.get_string(option_section, 'data_set', default='En n').split()
    ## Read section names for theory values
    config['theory_section']    = options.get_string(option_section, 'theory_section', default='cosebis')
    ## And the output name for the likelihood
    config['like_name'] = options.get_string(option_section, "like_name", default='cosebis')
    
    ## Load data & cov file
    try:
        TP_data= twopoint.TwoPointFile.from_fits(config['data_file'], covmat_name='COVMAT')
    except:
        raise OSError('\"%s\" not found' % config['data_file'])

    ## Extract the data vector & covariance matrix & put in config dict
    TP_data.choose_data_sets([config['stat_name']])
    config['data']       = TP_data.spectra[0].value
    config['covariance'] = TP_data.covmat
    config['inv_covariance'] = np.linalg.inv(TP_data.covmat)
        
    return config

def execute(block, config):

    ## Check that the theory has been calculated
    ## For cosebis this requires the cl_to_cosebis_interface to be run
    section_name=config['theory_section']
    if not block.has_section(section_name):
            raise AssertionError('I tried to grab theory values from \"%s\" but could not find anything.' % (section_name))
    angle = block[section_name, config['angle_name']]

    ## Building the theory vector to macth the ordering of the data vector
    sBuilder = SpectrumBuilder()
    nbTomo_max = 1000
    for i in range(nbTomo_max):
        if not block.has_value(section_name, 'bin_%d_%d' % (i+1, i+1)):
            break
                
        for j in range(i, nbTomo_max):
            value_name  = 'bin_%d_%d' % (j+1, i+1)
            try:
                value = block[section_name, value_name]
                sBuilder.addTomo(i, j, angle, value)
            except:
                break

    d = config['data']
    mu = np.concatenate(sBuilder.valList)   #Theory

    inv_cov = config['inv_covariance']
    r = d - mu

    chi2 = float(r @ inv_cov @ r)
    ln_like = -0.5*chi2

    like_name=config['like_name']
    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def clean(config):
    pass
