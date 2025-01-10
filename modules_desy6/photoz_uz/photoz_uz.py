from builtins import range
from cosmosis.datablock import option_section, names
import numpy as np

def setup(options):
    sample = options.get_string(option_section, "sample", "")
    pz = 'nz_'+config['sample']
    bias_section = sample+"_photoz_u"
    basis_file = options.get_string(option_section, "basis_file", "")
    n_modes = options.get_int(option_section, "n_modes", 0)
    z = block[pz, "z"]
    U = np.loadtxt(basis_file)[:,:n_modes].reshape((n_modes,block[pz, "nbin"],len(z)))
    perbin = options.get_string(option_section, "perbin", False)
    rescale = options.get_string(option_section, "rescale", False)
    return {"sample": sample,
            "bias_section": bias_section,
            "basis":U,
            "n_modes":n_modes,
            "rescale":rescale,
            "perbin":perbin}

def execute(block, config):
    pz = 'nz_'+config['sample']
    U = config['basis']
    n_modes = config['n_modes']
    nbin = block[pz, "nbin"]
    z = block[pz, "z"]
    uvals = config['bias_section']
    u_ = np.zeros((nbins,n_modes))
    for i in range(1,nbin+1):
        if perbin:
            u_[i-1,:] = np.array([block[uvals, "u_{0}_{1}".format(i-1,j) ] for j in range(n_modes)])
        else:
            u_[i-1,:] = np.array([block[uvals, "u_{0}".format(j) ] for j in range(n_modes)])
    if rescale:
        # replace below line with call to rescale the values from an initial unit gaussian.
        # Otherwise it will use the gaussian prior provided in the cosmosis file.
        # u_ is an array of the values for each bin with shape (nbins,n_modes).
        # if perbin=False, then the nbins dimension is just repeated entries.
        u_=u_
    for i in range(1,nbin+1):
        dz = U[:,i-1,:]@u_[i-1,:]
        bin_name = "bin_%d" % i
        nz = block[pz, bin_name]
        nz[1:]+=dz
        nz /= np.trapz(nz, z)
        block[pz, bin_name] = nz
    return 0

def cleanup(config):
    pass
