#coding:utf-8
#Module to compute growth factor from linear matter power spectrum
import numpy as np
from cosmosis.datablock import names, option_section
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

growth_section = "GROWTH_PARAMETERS"

def setup(options):
    k_growth = options.get_double(option_section, "k_growth", 1.e-3)
    return k_growth

def execute(block, config):

    k_growth=config
    z_pk,k_pk,p_lin = block.get_grid(names.matter_power_lin, 'z', 'k_h', 'p_k')
    a_pk = 1/(1+z_pk)
    growth_ind=np.where(k_pk>k_growth)[0][0]
    growth_array = np.sqrt(np.divide(p_lin[:,growth_ind], p_lin[0,growth_ind], 
                    out=np.zeros_like(p_lin[:,growth_ind]), where=p_lin[:,growth_ind]!=0.))
    #compute f(a) = dlnD/dlna, need to reverse when splining as expects increasing x
    growth_spline = IUS( np.log(a_pk[::-1]), np.log(growth_array[::-1]) )
    f_array = growth_spline.derivative()(np.log(a_pk[::-1]))[::-1]
    block[growth_section, 'd_z'] = growth_array
    block[growth_section, 'f_z'] = f_array
    block[growth_section, 'z'] = z_pk
    return 0








