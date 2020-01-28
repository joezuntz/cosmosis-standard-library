from cosmosis.datablock import names, option_section
import numpy as np
import os
import ctypes

# Names used
sigma_block_name = "sigma_r"

def set_vector(options, vmin, vmax, dv, vec):
    """Read a vector-valued parameter from the parameter file either directly or via min,max,n"""

    if options.has_value(option_section, vec):
        return np.array(options[option_section, vec])
    else:
        xmin = options[option_section, vmin]
        xmax = options[option_section, vmax]
        dx = options[option_section, dv]

        return np.arange(xmin, xmax, dx)

def setup(options):

    # We use a class as a namespace here to return config information

    class config:

        # Name of the section to get matter power from
        matter_power  = options.get_string(option_section, "matter_power", "matter_power_lin")

        # Ranges over which to calculate results
        z_vec = set_vector(options, "zmin", "zmax", "dz", "z")
        r_vec = None
        m_vec = None

        use_m = options.get_bool(option_section, "use_m", False)
        if use_m:
            m_vec = set_vector(options, "logmmin", "logmmax", "dlogm", "logm")
            print('using mass')
        else:
            r_vec = set_vector(options, "rmin", "rmax", "dr", "r")
            print('using radius')

        # compute dsigma2dlnm
        calc_dsig = options.get_bool(option_section, "calc_dsig", False)

        # Constants needed later
        rho_c = 2.775e11 #  rho_c/h^2
        rho_c_4pi_3 = rho_c * 4 * np.pi / 3.
        log_rho_c_4pi_3 = np.log10(rho_c_4pi_3)

        # printout settings
        verbose = int(options.get_bool(option_section, "verbose", False))

        print('m_vec = ', m_vec)
        print('r_vec = ', r_vec)
        print('z_vec = ', z_vec)

        ##############################################################################
        ########################## C WRAPPING ########################################
        ##############################################################################
        
        # C types
        C_VEC_DOUB = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        C_VEC_INT  = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
        C_INT      = ctypes.c_int
        C_DOUB     = ctypes.c_double

        # Importing these manually seems to be necessary on some systems
        # I thought it should be fine to do the linking in the Makefile
        # but apparently there is something wrong with that as written.
        # Would be nice to figure that out.
        dll1 = ctypes.CDLL('libgslcblas.so', mode=ctypes.RTLD_GLOBAL)
        dll2 = ctypes.CDLL('libgsl.so', mode=ctypes.RTLD_GLOBAL)

        # Import as Shared Object
        dirname = os.path.dirname(os.path.realpath(__file__))
        c_code = '{}/sigma.so'.format(dirname)
        lib_Sigma = ctypes.CDLL(c_code, mode=ctypes.RTLD_GLOBAL)

        # Specify the return and argument types of the function we will use
        Sigma_Func = lib_Sigma.executemain
        Sigma_Func.restype = C_INT
        Sigma_Func.argtypes = [
            C_DOUB        ,   # double* init_parameters
            C_VEC_INT     ,   # int*    int_config
            C_VEC_DOUB    ,   # double* Pk_k
            C_VEC_DOUB    ,   # double* Pk_z
            C_VEC_DOUB    ,   # double* Pk
            C_VEC_DOUB    ,   # double* z_vec
            C_VEC_DOUB    ,   # double* m_vec
            C_VEC_DOUB    ,   # double* r_vec
            C_VEC_DOUB    ,   # double* sigma2_m
            C_VEC_DOUB    ,   # double* dsigma2dlnM_m         
            ]


    return config

def execute(block, config):

    OmM = block[names.cosmological_parameters, "omega_m"]

    # Just a simple rename for clarity.
    z_vec = config.z_vec
    m_vec = config.m_vec
    r_vec = config.r_vec
    matter_power = config.matter_power

    if m_vec is None:
        m_vec = config.log_rho_c_4pi_3 + np.log10(OmM) + 3*np.log10(r_vec)

    if r_vec is None:
        r_vec = 0 * m_vec

    if config.verbose:
        print('m_vec', m_vec)
        print('z_vec', z_vec)

    nm, nz = len(m_vec), len(z_vec)
    
    # Load the matter power spectrum from the block, e.g. from CAMB.
    zk, k, Pk = block.get_grid(matter_power, "z", "k_h", "p_k")
    Pk = Pk.flatten()


    # For now OpenMP is disabled
    proc_num = 1

    # integer configuration parameters for the C++ code
    int_config = np.array([proc_num, len(k), nm, nz, len(zk), config.verbose, config.calc_dsig], dtype=np.int32)

    # Space where the output results will be stored
    sigma2_m = np.zeros(nm*nz)
    dsigma2dlnM_m = np.zeros(nm*nz)

    if config.verbose:
        print("************************** cluster code begin **********************************")
        print('len(Pk):',len(Pk))

    # Run the main function
    config.Sigma_Func (OmM, int_config, k, zk, Pk, z_vec, m_vec, r_vec, sigma2_m, dsigma2dlnM_m)

    if config.verbose:
        print("************************** cluster code ok **********************************")

    sigma2_m = sigma2_m.reshape((nz, nm))

    # We save the grid sigma2(R,z) and the logm vector separately.
    # This ordering is consistent with the old sigmar module
    block.put_grid(sigma_block_name, "R", r_vec, "z", z_vec, "sigma2", sigma2_m.T)
    block[sigma_block_name, "logm" ] = m_vec

    if config.calc_dsig:
        block[sigma_block_name, "dsigma2dlnM" ] = dsigma2dlnM_m.reshape(nz, nm).T

    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass



