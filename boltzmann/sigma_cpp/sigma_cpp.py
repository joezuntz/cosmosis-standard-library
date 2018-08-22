from cosmosis.datablock import names, option_section
import numpy as np
import os

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters
sigma = "sigma_r" #names.sigma_cpp

def set_vector(options, vmin, vmax, dv, vec):

    if options.has_value(option_section, vec):

        return np.array(options[option_section, vec])

    else:

        xmin = options[option_section, vmin]
        xmax = options[option_section, vmax]
        dx = options[option_section, dv]

        return np.arange(xmin, xmax, dx)

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.

    #Use this syntax to get a single parameter from the ini file section
    #for this module.  There is no type checking here - you get whatever the user
    #put in.


    #if options.has_value(option_section, "r"):

    class res():

        matter_power  = options[option_section, "matter_power"]

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

        rho_c = 2.775e11 #  rho_c/h^2
        rho_c_4pi_3 = rho_c * 4 * np.pi / 3.
        log_rho_c_4pi_3 = np.log10(rho_c_4pi_3)

        prt_detail = int(options.get_bool(
                            option_section, "prt_detail", False))

        print('m_vec = ', m_vec)
        print('r_vec = ', r_vec)
        print('z_vec = ', z_vec)

        ##############################################################################
        ########################## C WRAPPING ########################################
        ##############################################################################
        import ctypes
        
        # C types
        C_VEC_DOUB = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        C_VEC_INT  = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
        C_INT      = ctypes.c_int
        C_DOUB     = ctypes.c_double
        
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
            C_VEC_DOUB        # double* sigma_m
            ]

    #Now you have the input options you can do any useful preparation
    #you want.  Maybe load some data, or do a one-off calculation.

    #Whatever you return here will be saved by the system and the function below
    #will get it back.  You could return 0 if you won't need anything.


    return res

def execute(block, config):

    OmM = block[cosmo, "omega_m"]

    # Just a simple rename for clarity.
    z_vec = config.z_vec
    m_vec = config.m_vec
    r_vec = config.r_vec
    matter_power = config.matter_power

    if m_vec is None:
        m_vec = config.log_rho_c_4pi_3 + np.log10(OmM) + 3*np.log10(r_vec)

    if r_vec is None:
        r_vec = 0 * m_vec

    if config.prt_detail:
        print('m_vec', m_vec)
        print('z_vec', z_vec)
    
    zk  = block[matter_power, 'z']
    k  = block[matter_power, 'k_h']
    Pk = block[matter_power, 'p_k'].flatten()

    nm, nz = len(m_vec), len(z_vec)
    proc_num = 1

    int_config = np.array([proc_num, len(k), nm, nz, len(zk), config.prt_detail], dtype=np.int32)

    sigma_m = np.zeros(nm*nz)

    if config.prt_detail:
        print("************************** cluster code begin **********************************")
        print('len(Pk):',len(Pk))
    config.Sigma_Func (  
        OmM                 ,
        int_config          ,
        k                   ,
        zk                  ,
        Pk                  ,
        z_vec               ,
        m_vec               ,
        r_vec               ,
        sigma_m         
        )
    if config.prt_detail:
        print("************************** cluster code ok **********************************")

    sigma_m = sigma_m.reshape((nz, nm))

    # We save the grid sigma2(R,z) and the m vector separately.
    # This ordering is consistent with the old sigmar module
    block.put_grid(sigma, "R", r_vec, "z", z_vec, "sigma2", sigma_m.T)
    block[sigma, "m" ] = m_vec

    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass



