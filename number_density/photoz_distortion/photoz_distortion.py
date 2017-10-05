from builtins import range
from cosmosis.datablock import option_section, names
from scipy.interpolate import interp1d
import numpy as np

MODES = ["skew", "mean", "width"]


def setup(options):
    additive_bias = options[option_section, "mean"]
    broadening = options[option_section, "width"]
    skew = options[option_section, "tail"]
    per_bin = options[option_section, "bias_per_bin"]
    sample = options.get_string(option_section, "sample", default=None)
    if (not additive_bias) and (not broadening) and (not skew):
        raise ValueError("please set one or more of: %r to T" % MODES)

    try:
        cat_opt = options.get_string(option_section, "catastrophic_outliers")
    except:
        cat_opt = None

    return {"additive": additive_bias, "broadening": broadening, "skew": skew, "sample": sample, "per_bin": per_bin, "catastrophic_outliers": cat_opt}


def execute(block, config):
    additive, broadening, skew = config['additive'], config['broadening'], config['skew']
    pz = config['sample']
    biases = pz
    nbin = block[pz, "nzbin"]
    z = block[pz, "z"]

    # Pad the nz with zeros to prevent an unphysical cutoff
    # if the distribution is shifted upwards in redshift
    z_med = z[int(len(z) / 2.)]
    d = z[1] - z[0]
    add_pts = int((z[-1] - z_med) / d)
    pad = np.zeros(add_pts)
    padz = np.arange(z.max() + d, z.max() + (add_pts + 1) * d, d)
    z = np.append(z, padz)
    for i in range(1, nbin + 1):
        bin_name = "bin_%d" % i
        nz = block[pz, bin_name]
        nz = np.append(nz, pad)
        if config['per_bin']:
            if additive:
                bias = block[biases, "bias_%d" % i]
            if broadening:
                S = block[biases, "S_z_%d" % i]
            if skew:
                T = block[biases, "T_z_%d" % i]
        else:
            if additive:
                bias = block[biases, "bias_1"]
            if broadening:
                S_z = block[biases, "S_z_1"]
            if skew:
                T = block[biases, "T_z_1"]
        dz = np.zeros_like(z)
        f = interp1d(z, nz, kind='cubic', fill_value=0.0, bounds_error=False)
        # if mode=="multiplicative":
        #	nz_biased = f(z*(1-bias))
        if broadening:
            # Use the main peak of n(z) as a pivot point about which to distort n(z)
            zp = z[np.argwhere(nz == nz.max())[0][0]]
            dz += S * (z - zp)
        if additive:
            dz -= bias
        nz_biased = f(z + dz)

        if skew:
            # Redistribute proability upwards in redshift witout altering
            # the peak of the n(z)
            delta1 = T * (z - zp) + nz.max()
            delta2 = -1. * T * (z - zp) + nz.max()
            nz_biased += delta1 - delta2
            np.putmask(nz_biased, nz_biased < 0., 0.)

        # normalise
        nz_biased /= np.trapz(nz_biased, z)

        # Add a population of catastrophic outliers
        # See Hearin et al (2010)

        if config['catastrophic_outliers'] != None:
            cat_mode = block[pz, 'method']
            fcat = block[pz, 'fcat']  # 0.05
            dzcat = block[pz, 'dzcat']  # 0.129
            zcat0 = block[pz, 'zcat0']  # 0.65
            zcat = block[pz, 'zcat']  # 0.5
            sigcat = block[pz, 'sigcat']  # 0.1

            step = (dzcat / 2.) - abs(z - zcat0)
            step[step == 0.0] = -1.0
            step = 0.5 * (step / abs(step) + 1.0)

            dz = (z[1] - z[0])

            # Define a Gaussian island of outliers, normalised to
            # the probability scattered from the affected region
            if cat_mode == 'island':
                pcat = (1. / (2.0 * np.pi)**0.5 / sigcat) * np.exp(-1.0 *
                                                                   (z - zcat) * (z - zcat) / (2. * sigcat * sigcat))
                pcat *= np.trapz(step * fcat * nz_biased, z)
                nz_biased = (1. - step * fcat) * nz_biased + pcat

            # Or scatter it uniformly across the theory redshift range
            elif cat_mode == 'uniform':
                nz_biased = (1. - step * fcat) * nz_biased + \
                    step * fcat / (z[-1] - z[0])

        # renormalise
        nz_biased /= np.trapz(nz_biased, z)
        block[pz, bin_name] = nz_biased

    block[pz, 'z'] = z

    return 0


def cleanup(config):
    pass
