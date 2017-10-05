from __future__ import print_function
from numpy import log, pi, interp, where, loadtxt, dot
import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
likes = section_names.likelihoods
growthparams = section_names.growth_parameters
dist = section_names.distances

fsig_MEAN = 0.428  # Chuang et al 2013 BOSS DR9
fsig_SIGMA = 0.066
REDSHIFT = 0.57

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]
# Chuang et al 2013 BOSS DR9: H,D_a,omegamh2,bsigam8,fsigma8
DATA_file = os.path.join(ROOT_dir, 'chuang_etal_cmass_results.txt')
COV_file = os.path.join(ROOT_dir, 'chuang_etal_cmass_covariance.txt')

c_km_per_s = 299792.458


def setup(options):
    section = option_section
    mode = options.get_int(section, "mode", default=0)
    feedback = options.get_int(section, "feedback", default=0)
    if not mode:
        mean = options.get_double(section, "mean", default=fsig_MEAN)
        sigma = options.get_double(section, "sigma", default=fsig_SIGMA)
        redshift = options.get_double(section, "redshift", default=REDSHIFT)
        norm = 0.5 * log(2 * pi * sigma**2)
        return (mode, mean, sigma, norm, redshift, feedback)
    else:
        data = loadtxt(DATA_file)
        cov = loadtxt(COV_file)
        return (mode, data, cov, REDSHIFT, feedback)


def execute(block, config):
    # Configuration data, read from ini file above
    mode = config[0]
    if not mode:
        mean, sigma, norm, redshift, feedback = config[1:]
    else:
        data, cov, redshift, feedback = config[1:]

    z = block[growthparams, 'z']
    d_z = block[growthparams, 'd_z']
    f_z = block[growthparams, 'f_z']
    try:
        z0 = where(z == 0)[0][0]
    except IndexError:
        raise ValueError(
            "You need to calculate f(z) and d(z) down to z=0 to use the BOSS f*sigma8 likelihood")
    sig = block[cosmo, 'sigma_8']
    fsigma = (sig * (d_z / d_z[z0])) * f_z
    fsig = interp(redshift, z, fsigma)

    if feedback:
        print("Growth parameters: z = ", redshift, "fsigma_8  = ", fsig, " z0 = ", z0)

    if not mode:
        # compute the likelihood - just a simple Gaussian
        like = -(fsig - mean)**2 / sigma**2 / 2.0 - norm
        block[likes, 'BOSS_LIKE'] = like
        # signal that everything went fine
        return 0
    else:
        # distance info from datablock
        dist_z = block[dist, 'z']
        d_a = block[dist, 'd_a']
        h = c_km_per_s * block[dist, 'h']
        if dist_z[1] < dist_z[0]:
            dist_z = dist_z[::-1]
            d_a = d_a[::-1]
            h = h[::-1]

        h0 = block[cosmo, 'h0']
        omegam = block[cosmo, 'omega_m']
        bias = block[cosmo, 'bias']

        # sigma8 at z=0.57
        sigma_z = sig * (d_z / d_z[z0])
        s = interp(redshift, z, sigma_z)
        # D_a and H at z=0.57
        da_z = interp(redshift, dist_z, d_a)
        h_z = interp(redshift, dist_z, h)

        params = [h_z, da_z, omegam * h0**2, bias * s, fsig]
        d = params - data
        chi2 = dot(d, dot(cov, d))

        if feedback:
            print("[H(z=0.57),D_a(z=0.57),Omegamh2,b*sigma(z=0.57),fsigma8(z=0.57)] = [%lf,%lf,%lf,%lf,%lf]" %
                  (h_z, da_z, omegam * h0**2, bias * s, fsig))
        block[likes, 'BOSS_LIKE'] = chi2 * (-0.5)
        # signal that everything went fine
        return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness.  The joy of python.
    return 0
