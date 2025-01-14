from builtins import range
from cosmosis.datablock import option_section, names
import numpy as np

class Normalizer_cosmosis:
    def __init__(self, percentile, gmax, dg, xmax=10, kind='linear'):
        '''Class to create transformations from the gaussian 
        to the non-gaussian u distribution
        input: the percentile array of the initial u array
        generated with the same gmax, dg convention, 
        The methods  `degauss` will do  g->u.
        '''
        # Establishing matching percentile points
        gg = np.arange(-gmax,gmax+dg/2, dg)
        xx = np.copy(percentile)
        # Clean out any duplicates in the tables.  This is hacky.
        dup = xx[1:]==xx[:-1]
        keep = np.concatenate( (np.logical_not(dup), [True]))
        gg = gg[keep]
        xx = xx[keep]
        self.inv = interp1d(gg,xx,kind=kind, bounds_error=False, fill_value=(-xmax,xmax))
        return
    def degauss(self,g):
        '''Return de-Gaussianized values of g'''
        return self.inv(g)


def setup(options):
    sample = options.get_string(option_section, "sample", "")
    basis_file = options.get_string(option_section, "basis_file", "")
    #percentile_file = options.get_string(option_section, "Percentiles", "")
    n_modes = options.get_int(option_section, "n_modes", 0)
    n_bins = options.get_int(option_section, "n_bins", 0)
    #U = np.loadtxt(basis_file)[:,:n_modes]
    npzfile = np.load(basis_file)
    U=npzfile['U'][:,:n_modes]
    Percentiles=npzfile['Percentiles'][:,:n_modes]
    [gmax,dg,Ng]=npzfile['g_convention']
    lenz = len(U)//n_bins
    U = U.reshape((n_modes,n_bins,lenz))
    Percentiles = Percentiles.reshape((n_modes,n_bins,int(Ng)))
    #Percentiles = np.loadtxt(percentile_file)[:,:n_modes].reshape((n_modes,n_bins,Ng))
    Percentiles_classes=[[Normalizer_cosmosis(Percentiles[jmode][ibin], gmax, dg) for jmode in range(n_modes)] for ibin in range(n_bins)]
    perbin = options.get_string(option_section, "perbin", False)
    rescale = options.get_string(option_section, "rescale", False)
    return {"sample": sample,
            "bias_section": bias_section,
            "basis":U,
            "Percentiles_classes":Percentiles_classes,
            "n_modes":n_modes,
            "rescale":rescale,
            "perbin":perbin}

def execute(block, config):
    pz = 'nz_'+config['sample']
    U = config['basis']
    Percentiles_classes = config['Percentiles_classes']
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
        #u_=u_
        for ibin in range(nbin):
            for jmode in range(n_modes):
                u_[ibin,jmode]=  Percentiles_classes[ibin][jmode].degauss(u_[ibin,jmode])
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
