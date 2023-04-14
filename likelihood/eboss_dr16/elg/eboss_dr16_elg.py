import os
from numpy import log, pi, interp, where, loadtxt, dot, append, linalg, inf
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from scipy.interpolate import RegularGridInterpolator 
from scipy import zeros
from scipy import nan as scipynan

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
dist = section_names.distances

c_km_per_s =  299792.458

default_zeff = 0.845
default_rd_fiducial = 147.8

dirname = os.path.dirname(os.path.realpath(__file__))

def setup(options):
    
    section = option_section
    mode = options.get_int(section, "mode", default=False)
    feedback = options.get_bool(section, "feedback", default=False)
    
    zeff = options.get_double(section, "zeff", default_zeff)
    rd_fiducial = options.get_double(section, "rd_fiducial", default_rd_fiducial)

    # Read the relative probability table
    if not mode:
        print("ELG data")
        print("BAO only: Dv(z)/rd")
        data_file = os.path.join(dirname, 'sdss_DR16_ELG_BAO_DVtable.txt')
        data = loadtxt(data_file)
        out = [mode, feedback, zeff, rd_fiducial, data]
        
    else:
        print("ELG data")
        print("BAO+FS: Dm(z)/rd, Dh(z)/rd, and f(z)sigma8(z)")
        data_file = os.path.join(dirname, 'sdss_DR16_ELG_FSBAO_DMDHfs8gridlikelihood.txt')

        #from Arnaud de Mattia
        DM_npoints = 100
        DH_npoints = 100
        fs8_npoints = 100
        DM_grid = zeros(DM_npoints)
        DH_grid = zeros(DH_npoints)
        fs8_grid = zeros(fs8_npoints)
        prob_grid = zeros((DM_npoints,DH_npoints,fs8_npoints))
        tmp = loadtxt(data_file) 
        num = 0
        for i in range(DM_npoints):
            for j in range(DH_npoints):
                for k in range(fs8_npoints):
                    DM_grid[i],DH_grid[j],fs8_grid[k],prob_grid[i,j,k] = tmp[num]
                    num += 1

        loglike_spline = RegularGridInterpolator((DM_grid, DH_grid, fs8_grid), log(prob_grid), method='linear', bounds_error=False, fill_value=-inf)
        out = [mode, feedback, zeff, rd_fiducial, loglike_spline]     

    return out

def execute(block, config):
        
    # Redshift array
    z = block[dist, 'z']
    # Sound horizon at the drag epoch
    rd = block[dist, "rs_zdrag"]    
    # Comoving distance
    Dm = block[dist, 'd_m'] # in Mpc    
    # Hubble distance
    H = c_km_per_s*block[dist, 'H'] # in c/Mpc
    Dh = c_km_per_s/H
    # Spherically averaged distance    
    Dv = (z * Dm**2 *Dh)**(1/3) 
    
    # Configuration inputs    
    mode, feedback, zeff, rd_fiducial = config[0], config[1], config[2], config[3]
    
    # Find theory Dm, Dh and Dv at effective redshift by interpolation
    Dm_z_rd = interp(zeff, z, Dm)/rd
    Dh_z_rd = interp(zeff, z, Dh)/rd
    Dv_z_rd = interp(zeff, z, Dv)/rd
    
    # BAO-only
    if not mode:
        data = config[4]
        Dv_z_rd_data, like_data = data.T
        loglike = log(interp(Dv_z_rd, Dv_z_rd_data, like_data))
        if feedback:
            print()
            print('             zeff   pred')
            print('Dv_over_rd: %.3f  %.3f' % (zeff, Dv_z_rd))
            print('loglike:           %.3f' % (loglike))
            print()
    # BAO+FS
    else:
        z = block['growth_parameters', 'z']
        fsigma8 = block['growth_parameters', 'fsigma_8']
        # Find theory fsigma8 at fiducial redshift
        fsigma8_z = interp(zeff, z, fsigma8)
        
        loglike_spline = config[4]        
        loglike = loglike_spline([Dm_z_rd, Dh_z_rd, fsigma8_z])
        loglike = loglike[0]
        
        if feedback:
            print()
            print('             zeff   pred')
            print('Dm_over_rd: %.3f  %.3f' % (zeff, Dm_z_rd))
            print('Dh_over_rd: %.3f  %.3f' % (zeff, Dh_z_rd))
            print('loglike:           %.3f' % (loglike))
            print()

    block[likes,'eboss16_elg_like'] = loglike

    return 0

def cleanup(config):
    return 0

