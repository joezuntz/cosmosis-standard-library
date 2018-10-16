from cosmosis.datablock import option_section, names
import ctypes as ct
import numpy.ctypeslib as ctl
import numpy as np

import os
import sys

dirname = os.path.split(__file__)[0]
tp_dir = os.path.abspath(os.path.join(dirname,'../../likelihood/2pt/'))
sys.path.append(tp_dir)
import twopoint


# This opens a library written in C
lib = ct.cdll.LoadLibrary("{}/tessore_wigner.so".format(dirname))

# This loads a function from the C code and tells python its arguments
# and return value
lib.transform_cl_to_corr.restype = ct.c_int
lib.transform_cl_to_corr.argtypes = [
    ct.c_int,  # calc_type
    ct.c_int,  # ell_max
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags="C_CONTIGUOUS"),  # c_ell
    ct.c_int,  # n_theta
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags="C_CONTIGUOUS"),  #theta
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C_CONTIGUOUS","WRITEABLE"))  #xi
]

# These are constants that are defined in the C code and passed to it
# to describe which transform to do
XI_PLUS = ct.c_int.in_dll(lib, "XI_PLUS").value
XI_MINUS = ct.c_int.in_dll(lib, "XI_MINUS").value
GAMMA_T = ct.c_int.in_dll(lib, "GAMMA_T").value
W_THETA = ct.c_int.in_dll(lib, "W_THETA").value


XIPLUS_TYPE  = (twopoint.Types.galaxy_shear_plus_real, twopoint.Types.galaxy_shear_plus_real)
XIMINUS_TYPE = (twopoint.Types.galaxy_shear_minus_real, twopoint.Types.galaxy_shear_minus_real)
GAMMAT_TYPE  = (twopoint.Types.galaxy_position_real, twopoint.Types.galaxy_shear_plus_real)
WTHETA_TYPE  = (twopoint.Types.galaxy_position_real, twopoint.Types.galaxy_position_real)


def transform(calc_type, ell, c_ell, theta):
    # Check calculation type
    if calc_type not in [XI_PLUS, XI_MINUS, GAMMA_T, W_THETA]:
        raise ValueError("Parameter calc_type must be 0 (XI_PLUS), 1 (XI_MINUS), 2 (GAMMA_T), 3 (W_THETA)")

    # We need every ell value from 0 to ell_max, so interpolate to get them
    ell_max = int(ell.max())
    interpolated_ell = np.arange(ell_max+1, dtype=int)
    interpolated_cell = np.interp(interpolated_ell, ell, c_ell, left=0., right=0.)

    # Make sure theta is an array
    theta = np.array(theta, dtype=np.double)
    ntheta = len(theta)

    # Set up output result space
    xi = np.zeros(ntheta, dtype=np.double)

    # Call C transform code.
    status = lib.transform_cl_to_corr(calc_type, ell_max, interpolated_cell, ntheta, theta, xi)

    return xi

def setup_from_file(options):
    raise NotImplementedError("setup_from_file not finished. email joe if you want to help out.")
    
    filename = options[option_section, "twopoint_file"]
    data_names = options[option_section, "data_names"].split()
    T=twopoint.TwoPointFile.from_fits(filename)

    corr_types = []
    input_sections = []
    output_sections = []
    theta_vals = []

    for data_name in data_names:
        S=T.get_spectrum(data_name)
        types = (S.type1, S.type2)
        try:
            corr_type = {
                XIPLUS_TYPE: XI_PLUS,
                XIMINUS_TYPE: XI_MINUS,
                GAMMAT_TYPE: GAMMA_T,
                WTHETA_TYPE: W_THETA
            }[types]
            corr_types.append(corr_type)
        except KeyError:
            raise ValueError("Error: could not understand type of data {} in file {}".format(data_name,filename))

        input_section = {
            XIPLUS_TYPE: "shear_cl",
            XIMINUS_TYPE: "shear_cl",
            GAMMAT_TYPE: "galaxy_shear_cl",
            WTHETA_TYPE: "galaxy_cl",
        }
        
        output_section = {
            XIPLUS_TYPE: "shear_xi",
            XIMINUS_TYPE: "shear_xi",
            GAMMAT_TYPE: "galaxy_shear_xi",
            WTHETA_TYPE: "galaxy_xi",
        }
        
        theta = arg_this_does_not_work


def setup_from_values(options):
    corr_type = options[option_section, "corr_type"]

    if corr_type==0 or corr_type=="xi":
        corr_types=[XI_PLUS, XI_MINUS]
        default_input = "shear_cl"
        default_output = "shear_xi"
        n = 2
    elif corr_type==1 or corr_type=="gammat":
        corr_types=[GAMMA_T]
        default_input = "galaxy_shear_cl"
        default_output = "galaxy_shear_xi"
        n = 1
    elif corr_type==2 or corr_type=="wtheta":
        corr_types=[W_THETA]
        default_input = "galaxy_cl"
        default_output = "galaxy_xi"
        n = 1
    else:
        raise ValueError("Unknown value of corr_type parameter. Should be 0/xi  1/gammat 2/wtheta")

    input_sections = [options.get_string(option_section, "input_section_name", default=default_input)]*n
    output_sections = [options.get_string(option_section, "output_section_name", default=default_output)]*n


    theta_min = options.get_double(option_section, "theta_min")
    theta_max = options.get_double(option_section, "theta_max")
    n_theta = options.get_int(option_section, "n_theta")

    # Make theta in radians
    theta_vals = [
        np.radians(np.linspace(theta_min, theta_max, n_theta)/60.),
    ]*n


    return (corr_types, input_sections, output_sections, theta_vals)

def setup(options):

    if options.has_value(option_section, "twopoint_file"):
        config = setup_from_file(options)
    else:
        config = setup_from_values(options)

    return config



def execute(block, config):

    (corr_types, input_sections, output_sections, thetas) = config

    n = len(corr_types)

    for i,(corr_type,input_section,output_section,theta) in enumerate(zip(*config)):
        if corr_type == XI_PLUS:
            out_name = "xiplus_{}_{}"
        elif corr_type == XI_MINUS:
            out_name = "ximinus_{}_{}"
        else:
            out_name = "bin_{}_{}"

        ell = block[input_section, "ell"]
        nbin_a = block[input_section, "nbin_a"]
        nbin_b = block[input_section, "nbin_b"]

        block[output_section, "theta"] = theta

        for b1 in range(1,nbin_a+1):
            for b2 in range(1,nbin_b+1):
                name = "bin_{}_{}".format(b1,b2)
                if not block.has_value(input_section, name):
                    continue
                c_ell = block[input_section, name]
                xi = transform(corr_type, ell, c_ell, theta)
                block[output_section, out_name.format(b1,b2)] = xi
                print("Saving {} {}".format(output_section, out_name.format(b1,b2)))
    return 0




def test():
    xis = []
    n_ell = 5000
    ell = np.logspace(0.,7.,n_ell)
    c_ell = 1./ell**2
    theta = np.logspace(0., 3., 50) / 60.
    xi_5000 = transform(W_THETA, ell, c_ell, theta)
    print("Done full")
    for ell_max in [200000, 400000, 600000, 800000, 1000000]:
        for n_ell in [50,100,200,300,400,500]:
            ell = np.logspace(0.,np.log10(ell_max),n_ell)
            c_ell = 1./ell**2
            xi = transform(W_THETA, ell, c_ell, theta)
            r = abs((xi/xi_5000)-1)
            rmax=r.max()
            print(f'nell={n_ell},lmax={ell_max}  rmax = {rmax:.2e}')


if __name__ == '__main__':
    test()