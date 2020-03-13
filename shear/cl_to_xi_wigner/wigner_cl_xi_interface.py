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
lib = ct.cdll.LoadLibrary("{}/wigner_d.so".format(dirname))

# This loads a function from the C code and tells python its arguments
# and return value
# void wigner_d(int l0, int l1, int n, int m, double theta, double* d);
lib.wigner_d.restype = None
lib.wigner_d.argtypes = [ ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_double,
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C_CONTIGUOUS","WRITEABLE")) ]

# These are constants that describe which transform to do
XI_PLUS  = 0
XI_MINUS = 1
GAMMA_T  = 2
W_THETA  = 3

XIPLUS_TYPE  = (twopoint.Types.galaxy_shear_plus_real, twopoint.Types.galaxy_shear_plus_real)
XIMINUS_TYPE = (twopoint.Types.galaxy_shear_minus_real, twopoint.Types.galaxy_shear_minus_real)
GAMMAT_TYPE  = (twopoint.Types.galaxy_position_real, twopoint.Types.galaxy_shear_plus_real)
WTHETA_TYPE  = (twopoint.Types.galaxy_position_real, twopoint.Types.galaxy_position_real)

def transform(calc_type, ell, theta):
    # Check calculation type
    if calc_type not in [XI_PLUS, XI_MINUS, GAMMA_T, W_THETA]:
        raise ValueError("Parameter calc_type must be 0 (XI_PLUS), 1 (XI_MINUS), 2 (GAMMA_T), 3 (W_THETA)")

    # spins of the transform
    s1, s2 = {
        XI_PLUS:  (2, 2),
        XI_MINUS: (2,-2),
        GAMMA_T:  (2, 0),
        W_THETA:  (0, 0)
    }[calc_type]

    # We need every integer ell value from lmin to lmax (inclusive)
    lmin = np.max([np.abs(s1), np.abs(s2)])
    lmax = int(np.max(ell))
    lint = np.arange(lmin, lmax+1)

    # Indices corresponding to ell values
    lidx = np.arange(len(ell))

    # transformation matrix to go from cl to xi
    tfm = np.zeros((len(theta), len(ell)))

    # Go through theta values
    for i, t in enumerate(theta):

        # Set up output result space
        d = np.zeros(lmax-lmin+1, dtype=np.double)

        # Call C code.
        lib.wigner_d(lmin, lmax, s1, s2, t, d)

        # Create the transformation from all C_l to xi of this theta:
        # Summing over all integer lmin <= l <= lmax, the terms would be
        #     (2*l + 1)/(4*pi)*C[l]*d[l]
        # However, C[l] at integer l is interpolated from the array of values
        # at non-integer ell:
        #     C[l] = (1-x)*c_ell[k] + x*c_ell[k+1]
        # Hence collect the contribution of each c_ell to xi in a matrix
        for l, z in zip(lint, (2*lint + 1)/(4*np.pi)*d):
            # position of l in ell array
            x = np.interp(l, ell, lidx, left=-1, right=-1)

            # integer and fractional part for interpolation
            k = int(x)
            x = x - k

            # check for hitting the last value
            if k == lidx[-1]:
                k = k-1
                x = 1

            # interpolate
            if k != -1:
                tfm[i,k] += (1-x)*z
                tfm[i,k+1] += x*z

    # the array to transform from c_ell to theta
    return tfm

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
    elif corr_type==1 or corr_type=="wtheta":
        corr_types=[W_THETA]
        default_input = "galaxy_cl"
        default_output = "galaxy_xi"
        n = 1
    elif corr_type==2 or corr_type=="gammat":
        corr_types=[GAMMA_T]
        default_input = "galaxy_shear_cl"
        default_output = "galaxy_shear_xi"
        n = 1
    else:
        raise ValueError("Unknown value of corr_type parameter. Should be 0/xi, 1/wtheta, 2/gammat")

    input_sections = [options.get_string(option_section, "input_section_name", default=default_input)]*n
    output_sections = [options.get_string(option_section, "output_section_name", default=default_output)]*n

    theta_min = options.get_double(option_section, "theta_min")
    theta_max = options.get_double(option_section, "theta_max")
    n_theta = options.get_int(option_section, "n_theta")

    # Make theta in radians
    theta_vals = [
        np.radians(np.linspace(theta_min, theta_max, n_theta)/60.),
    ]*n

    # Transformation matrices
    tfms = [[] for i in range(n)]

    return corr_types, input_sections, output_sections, theta_vals, tfms

def setup(options):
    if options.has_value(option_section, "twopoint_file"):
        config = setup_from_file(options)
    else:
        config = setup_from_values(options)

    return config

def execute(block, config):
    for i,(corr_type,input_section,output_section,theta,tfm) in enumerate(zip(*config)):
        if corr_type == XI_PLUS:
            output_section = output_section + '_plus'
        elif corr_type == XI_MINUS:
            output_section = output_section + '_minus'

        ell = block[input_section, "ell"]
        nbin_a = block[input_section, "nbin_a"]
        nbin_b = block[input_section, "nbin_b"]

        block[output_section, "theta"] = theta

        if len(tfm) == 0:
            print('computing transform type {} from ell ({:g}, {:g}, {}) to theta ({:g}, {:g}, {})'.format(corr_type, ell[0], ell[-1], len(ell), theta[0], theta[-1], len(theta)))
            tfm.append(transform(corr_type, ell, theta))

        for b1 in range(1,nbin_a+1):
            for b2 in range(1,nbin_b+1):
                name = "bin_{}_{}".format(b1,b2)
                if not block.has_value(input_section, name):
                    continue
                c_ell = block[input_section, name]
                xi = np.dot(tfm[0], c_ell)
                block[output_section, name] = xi
    return 0
