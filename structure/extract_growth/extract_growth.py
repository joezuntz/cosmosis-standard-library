from cosmosis.datablock import names, option_section
import scipy.interpolate
import numpy as np

def setup(options):
    input_section = options.get_string(option_section, "input_section", names.matter_power_lin)
    output_section = options.get_string(option_section, "output_section", names.growth_parameters)
    return [input_section, output_section]

def execute(block, config):
    input_section, output_section = config

    z,k,p = block.get_grid(input_section, "z", "k_h", "p_k")

    k0 = k.argmin()
    z0 = z.argmin()

    p_z = p[:, k0]

    D = np.sqrt(p_z / p_z[z0])
    logD = np.log(D)

    a = 1.0/(1.0+z)
    loga = np.log(a)

    logD_spline = scipy.interpolate.UnivariateSpline(loga[::-1], logD[::-1])
    f_spline = logD_spline.derivative()
    f = f_spline(loga)

    block[output_section, "z"] = z
    block[output_section, "d_z"] = D
    block[output_section, "f_z"] = f

    return 0