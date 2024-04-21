from builtins import range
from cosmosis.datablock import option_section, names
import numpy as np

def setup(options):
    sample = options.get_string(option_section, "sample", "")
    bias_section = sample+"_photoz_u"
    basis_file = options.get_string(option_section, "basis_file", "")
    n_modes = options.get_int(option_section, "n_modes", 0)
    U = np.loadtxt(basis_file)[:,:n_modes]
    return {"sample": sample,
            "bias_section": bias_section,
            "basis":U,
            "n_modes":n_modes}

def execute(block, config):
    pz = 'nz_'+config['sample']
    uvals = config['bias_section']
    U = config['basis']
    n_modes = config['n_modes']
    nbin = block[pz, "nbin"]
    if nbin>4:
        nbin=4
    z = block[pz, "z"]
    dz = U@np.array([block[uvals, "u_%d" % j] for j in range(n_modes)])
    dz = dz.reshape((nbin,len(dz)//nbin))
    for i in range(1,nbin+1):
        bin_name = "bin_%d" % i
        nz = block[pz, bin_name]
        nz[1:]+=dz[i-1,:]
        nz /= np.trapz(nz, z)
        block[pz, bin_name] = nz
    return 0

def cleanup(config):
    pass
