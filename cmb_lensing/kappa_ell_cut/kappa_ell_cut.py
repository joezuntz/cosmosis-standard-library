from cosmosis.datablock import option_section, names
import numpy as np
import warnings

def get_nbins(block, section):
    if block.has_value(section, "nbin_a"):
        n_a=block[section,"nbin_a"]
        n_b=block[section,"nbin_b"]
    else:
        n_a=block[section,"nbin"]
        n_b=n_a
    return n_a, n_b


def apply_lcut(block, section, Lmin, Lmax):
    n_other, n_cmb = get_nbins(block, section)
    if (n_cmb != 1):
        raise ValueError("Found more than 1 CMB bin - there is definitely only 1 CMB.")

    ell = block[section, 'ell']
    w   = np.ones(ell.shape[0])
    w[ell<Lmin]=0
    w[ell>Lmax+1]=0
    for bini in xrange(0,n_other):
        C_ell_orig = block[section, 'bin_{}_1'.format(bini+1)]
        beamed_C_ell = C_ell_orig*w
        block[section, 'bin_{}_1'.format(bini+1)] = beamed_C_ell

def setup(options):
    shearkappa_section = options.get_string(option_section, "shearkappa_section", default="shear_cmbkappa_cl")
    galkappa_section = options.get_string(option_section, "galkappa_section", default="galaxy_cmbkappa_cl")

    if options.has_value(option_section,'Lmin'):
        Lmin = options[option_section, 'Lmin']
    else:
        warnings.warn("No Lmin supplied! Setting Lmin=0.")
        Lmin = 0
    if options.has_value(option_section,'Lmax'):
        Lmax = options[option_section, 'Lmax']
    else:
        warnings.warn("No Lmax supplied! Setting Lmax=999999.")
        Lmax = 999999

    return {"shearkappa_section":shearkappa_section, "galkappa_section":galkappa_section, "Lmin":Lmin, "Lmax":Lmax}

def execute(block, config):
    shearkappa_section = config['shearkappa_section']
    galkappa_section = config['galkappa_section']
    Lmin = config['Lmin']
    Lmax = config['Lmax']

    if block.has_section(galkappa_section):
        apply_lcut(block, galkappa_section, Lmin, Lmax)


    if block.has_section(shearkappa_section):
        apply_lcut(block, shearkappa_section, Lmin, Lmax)
        
    return 0

def cleanup(config):
    pass
