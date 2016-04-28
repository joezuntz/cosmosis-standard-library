from astropy.io import fits
import astropy.units
from astropy.table import Table
from enum34 import Enum
import numpy as np

#FITS header keyword indicating 2-point data
TWOPOINT_SENTINEL = "2PTDATA"
NZ_SENTINEL = "NZDATA"
COV_SENTINEL = "COVDATA"
window_types=["SAMPLE","CLBP"]



#Please do not add things to this list
ANGULAR_UNIT_TYPES=[
    astropy.units.arcsec,
    astropy.units.arcmin,
    astropy.units.rad,
    astropy.units.deg,
]
ANGULAR_UNITS={unit.name:unit for unit in ANGULAR_UNIT_TYPES}



class Types(Enum):
    """
    This is an enumeration - a list of possible values with names and code values
    that we can use in FITS headers.

    It enumerates the different quantities that a two-point measurement can correlate.
    For example, CMB T,E,B, lensing E and B, galaxy position and magnification.

    It specifies the quantity and whether we are in Fourier space or real space.

    One special case is xi_{-} and xi_{+} in galaxy shear.  These values are already
    correlations of combinations of E and B modes. We denote this as xi_{++} and xi_{--} here.

    """
    galaxy_position_fourier = "GPF"
    galaxy_shear_emode_fourier = "GEF"
    galaxy_shear_bmode_fourier = "GBF"
    galaxy_position_real = "GPR"
    galaxy_shear_plus_real = "G+R"
    galaxy_shear_minus_real = "G-R"

    @classmethod
    def lookup(cls, value):
        for T in cls:
            if T.value==value:
                return T



class NumberDensity(object):
    """
    This class contains n(z) information for a collection of numbered redshift bins.
    It is expected to be used for a single sample type (e.g. galaxy sample) split into
    tomographic bins rather than the more complex cases where there are multiple different
    quantities (e.g. two galaxy samples).

    Since the main use case for this is photometric redshifts, and photo-z codes typically
    produce histogram type data sets (that is, they look like step functions between each bin),
    that is what this form assumes.
    """
    def __init__(self, name, zlow, z, zhigh, nzs):
        self.name=name
        self.zlow = zlow
        self.z = z
        self.zhigh = zhigh
        self.nbin = len(nzs)
        if self.nbin>0:
            self.nsample = len(nzs[0])
        else:
            self.nsample = 0
        self.nzs = nzs

    @classmethod
    def from_fits(cls, extension):
        #load in the n(z) data
        data = extension.data

        z = data['Z_MID']
        zlow = data['Z_LOW']
        zhigh = data['Z_HIGH']
        i = 1
        name = 'BIN{}'.format(i)
        nzs = []
        while name in data.names:
            nz = data[name]
            nzs.append(nz)
            i+=1 
            name = 'BIN{}'.format(i)

        N = cls(extension.name, zlow, z, zhigh, nzs)

        return N

    def to_fits(self):
        header = fits.Header()
        header[NZ_SENTINEL]= True
        header['EXTNAME'] = self.name

        columns = [
            fits.Column(name='Z_LOW', array=self.zlow, format='D'),
            fits.Column(name='Z_MID', array=self.z, format='D'),
            fits.Column(name='Z_HIGH', array=self.zhigh, format='D'),
        ]

        for i in xrange(self.nbin):
            name = "BIN{}".format(i+1)
            columns.append(fits.Column(name=name, array=self.nzs[i], format='D'))

        extension = fits.BinTableHDU.from_columns(columns, header=header)
        return extension        
        


class SpectrumMeasurement(object):
    def __init__(self, name, bins, types, kernels, windows, angular_bin, value, 
        angle=None, error=None, angle_unit=None, metadata=None):
        """metadata is a dictionary which will get added to the fits header"""
        self.name = name
        self.bin1, self.bin2 = bins
        self.bin_pairs = self.get_bin_pairs()  #unique bin pairs
        self.type1, self.type2 = types
        self.kernel1, self.kernel2 = kernels
        self.angular_bin = angular_bin
        self.angle = angle
        self.value = value
        if windows in window_types:
            self.windows = windows
        else:
            raise TypeError("window type %s not recognised"%windows)
        self.error = error
        self.metadata = metadata
        if self.is_real_space():
            #angle is real
            msg = "Files with real-space units must specify units as one of: {}".format(ANGULAR_UNITS.keys())
            assert angle_unit in ANGULAR_UNITS,  msg
        self.angle_unit = angle_unit

    def get_bin_pairs(self):
        all_pairs = zip(self.bin1,self.bin2)
        unique_pairs=[]
        for p in all_pairs:
            if p not in unique_pairs:
                unique_pairs.append(p)
        return unique_pairs

    def is_real_space(self):
        return self.type1.value.endswith("R") or self.type2.value.endswith("R")

    def convert_angular_units(self, unit):
        if not self.is_real_space():
            raise ValueError("Two point spectrum has no units to convert; it is in Fourier space")

        if self.windows not in ["SAMPLE", "CLBP"]:
            raise NotImplementedError("Need to write code to transform window function units")
        old_unit = ANGULAR_UNITS[self.angle_unit]
        new_unit = ANGULAR_UNITS[unit]

        print "Converting angle units of {} from {} -> {} (factor {})".format(
            self.name,old_unit, new_unit, old_unit.to(new_unit))
        angle_with_units = self.angle*old_unit
        self.angle = angle_with_units.to(new_unit).value
        self.angle_unit = unit


    def apply_mask(self, mask):
        """mask is a boolean array which True for elements to be kept"""
        self.bin1 = self.bin1[mask]
        self.bin2 = self.bin2[mask]
        self.angular_bin = self.angular_bin[mask]
        self.angle = self.angle[mask]
        self.value = self.value[mask]

    def auto_bins(self):
        return self.bin1==self.bin2


    def __len__(self):
        return len(self.value)

    def nbin(self):
        return np.max([self.bin1.max(), self.bin2.max()])

    def get_pair(self, bin1, bin2):
        w = (self.bin1==bin1) & (self.bin2==bin2)
        return self.angle[w], self.value[w]

    def get_pair_mask(self, bin1, bin2):
        w = (self.bin1==bin1) & (self.bin2==bin2)
        return w

    def get_error(self, bin1, bin2):
        if self.error is None:
            return None
        w = (self.bin1==bin1) & (self.bin2==bin2)
        return self.error[w]

    @classmethod
    def from_fits(cls, extension, covmat_info=None):
        name=extension.name
        #determine the type of the quantity involved
        type1 = Types.lookup(extension.header['QUANT1'])
        type2 = Types.lookup(extension.header['QUANT2'])

        #and the name of the kernels and the window function
        #extensions
        kernel1 = extension.header['KERNEL_1']
        kernel2 = extension.header['KERNEL_2']
        windows = extension.header['WINDOWS']

        if windows not in window_types:
            raise NotImplementedError("Have not yet coded window functions for angular bins")

        #Now load the data
        data = extension.data
        bin1 = data['BIN1']
        bin2 = data['BIN2']
        angular_bin = data['ANGBIN']
        value = data['VALUE']
        if "ANG" in data.names:
            angle = data['ANG']
            ang_index = data.names.index("ANG")
            angle_unit= extension.header.get('TUNIT{}'.format(ang_index+1))
            if angle_unit is not None:
                angle_unit = angle_unit.strip()
        else:
            angle = None
            angle_unit = None


        #Load a chunk of the covariance matrix too if present.
        if covmat_info is None:
            error = None
        else:
            error = covmat_info.get_error(name)

        return SpectrumMeasurement(name, (bin1, bin2), (type1, type2), (kernel1, kernel2), windows,
            angular_bin, value, angle, error, angle_unit=angle_unit)

    def to_fits(self):
        header = fits.Header()
        header[TWOPOINT_SENTINEL]= True
        header['EXTNAME']=self.name
        header['QUANT1'] = self.type1.value
        header['QUANT2'] = self.type2.value
        header['KERNEL_1'] = self.kernel1
        header['KERNEL_2'] = self.kernel2
        header['WINDOWS'] = self.windows #NOT YET CODED ANYTHING ELSE
        header['N_ZBIN_1'] = len(np.unique(self.bin1))
        header['N_ZBIN_2'] = len(np.unique(self.bin2))
        if self.metadata is not None:
            #Check metadata doesn't share any keys with the stuff that's already in the header
            assert set(self.metadata.keys()).isdisjoint(header.keys())
            for key,val in self.metadata.iteritems():
                header[key]=val
        header['N_ANG']=len(np.unique(self.angular_bin))

        columns = [
            fits.Column(name='BIN1', array=self.bin1, format='K'),
            fits.Column(name='BIN2', array=self.bin2, format='K'),
            fits.Column(name='ANGBIN', array=self.angular_bin, format='K'),
            fits.Column(name='VALUE', array=self.value, format='D'),
        ]
        if self.angle is not None:
            if self.windows=="SAMPLE":
                columns.append(fits.Column(name='ANG', array=self.angle, format='D', unit=self.angle_unit))
            if self.windows=="CLBP":
                columns.append(fits.Column(name='ANG', array=self.angle, format='2K',unit=self.angle_unit))

        extension = fits.BinTableHDU.from_columns(columns, header=header)
        return extension



class CovarianceMatrixInfo(object):
    """Encapsulate a covariance matrix and indices in it"""
    def __init__(self, name, names, lengths, covmat):
        super(CovarianceMatrixInfo, self).__init__()
        self.name = name
        self.names = names
        self.lengths = lengths
        self.starts=[0]
        for i,l in enumerate(self.lengths[:-1]):
            self.starts.append(l+self.starts[i])
        self.covmat = covmat
        self.diagonal = covmat.diagonal()

    def get_error(self, name):
        i = self.names.index(name)
        start = self.starts[i]
        end = start + self.lengths[i]
        return self.diagonal[start:end]**0.5

    def to_fits(self):
        header = fits.Header()
        header[COV_SENTINEL]= True
        header['EXTNAME'] = self.name
        for i,(start_index,name) in enumerate(zip(self.starts, self.names)):
            header['STRT_{}'.format(i)] = start_index
            header['NAME_{}'.format(i)] = name
        extension=fits.ImageHDU(data=self.covmat, header=header)
        return extension


    @classmethod
    def from_fits(cls, extension):
        cov_name = extension.name
        covmat = extension.data
        header = extension.header
        i = 0
        measurement_names = []
        start_indices = []
        while True:
            name_card = 'NAME_{}'.format(i)
            if name_card not in header:
                break
            measurement_names.append(header[name_card])
            start_index_card = 'STRT_{}'.format(i)
            start_indices.append(header[start_index_card])
            i+=1
        lengths = []
        current_length = 0
        #this only works if more than one spectrum
        if len(start_indices)>1:
            for start, end in zip(start_indices[:-1], start_indices[1:]):
                lengths.append(end-start)
            if start_indices:
                lengths.append(covmat.shape[0]-start_indices[-1])
        else:
            lengths.append(covmat.shape[0])
        return cls(cov_name, measurement_names, lengths, covmat)




class TwoPointFile(object):
    def __init__(self, spectra, kernels, windows, covmat_info):
        if windows is None:
            windows = {}
        self.spectra = spectra
        self.kernels = kernels
        self.windows = windows
        self.covmat_info = covmat_info
        if covmat_info:
            #self.covmat = covmat_info.covmat
            self.covmat = self.get_cov_start()
        else:
            self.covmat = None

    def get_spectrum(self, name):
        spectra = [spectrum for spectrum in self.spectra if spectrum.name==name]
        n = len(spectra)
        if n==0:
            raise ValueError("Spectrum with name %s not found in file"%name)
        elif n>1:
            raise ValueError("Multiple spectra with name %s found in file"%name)
        else:
            return spectra[0]

    def get_kernel(self, name):
        return self.kernels[name]

    def _mask_covmat(self, masks):
        #Also cut down the covariance matrix accordingly
        if self.covmat is not None:
            mask = np.concatenate(masks)
            self.covmat = self.covmat[mask, :][:, mask]

    def mask_bad(self, bad_value):
        "Go through all the spectra masking out data where they are equal to bad_value"
        masks = []
        #go through the spectra and covmat, masking out the bad values.
        for spectrum in self.spectra:
            #nb this will not work for NaN!
            mask = (spectrum.value != bad_value) 
            spectrum.apply_mask(mask)
            print "Masking {} values in {}".format(mask.size-mask.sum(), spectrum.name)
            #record the mask vector as we will need it to mask the covmat
            masks.append(mask)
        if masks:
            self._mask_covmat(masks)

    def mask_cross(self):
        masks = []
        for spectrum in self.spectra:
            mask = spectrum.auto_bins()
            spectrum.apply_mask(mask)
            print "Masking {} cross-values in {}".format(mask.size-mask.sum(), spectrum.name)
            masks.append(mask)
        if masks:
            self._mask_covmat(masks)

    def mask_scales(self, cuts={}, bin_cuts=[]):
        masks=[]
        print
        for spectrum in self.spectra:
            mask = np.ones(len(spectrum), dtype=bool)
            for b1,b2 in spectrum.bin_pairs:
                w_full = np.where((spectrum.bin1==b1) & (spectrum.bin2==b2))[0]
                if (spectrum.name, b1, b2) in bin_cuts:
                    print "Removing {} bin ({},{}) altogether.".format(spectrum.name, b1, b2)
                    mask[w_full] = False
                    continue

                cut = cuts.get((spectrum.name,b1,b2))
                if cut is None:
                    print "No cut specified for {} bin ({},{})".format(spectrum.name, b1, b2)
                    continue

                #Actually do the cut
                ang_min, ang_max = cut
                w = np.where((spectrum.bin1==b1) & (spectrum.bin2==b2) & 
                    ((spectrum.angle<ang_min) | (spectrum.angle>ang_max) ) )[0]
                
                print "Cutting {} bin pair ({},{}) to angle range ({} - {}) : this removes {} values out of {}".format(
                    spectrum.name, b1, b2, ang_min, ang_max, len(w), len(w_full))

                mask[w] = False
            masks.append(mask)
            spectrum.apply_mask(mask)
            print

        if masks:
            self._mask_covmat(masks)
            
            
    def mask_scale(self, spectra_to_cut, min_scale=-np.inf, max_scale=np.inf):
        masks = []
        #go through the spectra and covmat, masking out the bad values.
        for spectrum in self.spectra:
            if (spectra_to_cut!="all") and (spectrum.name not in spectra_to_cut):
                continue
            #nb this will not work for NaN!
            mask = (spectrum.angle > min_scale) & (spectrum.angle < max_scale) 
            spectrum.apply_mask(mask)            
            print "Masking {} values in {} because they had ell or theta outside ({},{})".format(mask.size-mask.sum(), spectrum.name, min_scale, max_scale)
            #record the mask vector as we will need it to mask the covmat
            masks.append(mask)

        if masks:
            self._mask_covmat(masks)

    def choose_data_sets(self, data_sets):
        """Strip out any data sets not in the given list."""
        data_sets = [d.lower() for d in data_sets]
        mask = []
        use = []
        for spectrum in self.spectra:
            if spectrum.name.lower() in data_sets:
                use.append(True)
                mask.append(np.ones(spectrum.bin1.size, dtype=bool))
            else:
                use.append(False)
                mask.append(np.zeros(spectrum.bin1.size, dtype=bool))
        for data_set in data_sets:
            if not any(spectrum.name.lower()==data_set for spectrum in self.spectra):
                raise ValueError("Data set called {} not found in two-point data file.".format(data_set))
        self.spectra = [s for (u,s) in zip(use,self.spectra) if u]
        if self.covmat is not None:
            mask = np.concatenate(mask)
            self.covmat = self.covmat[mask,:][:,mask]

    def get_cov_start(self):

        #This gets the covariance array in the right order (before any scale cuts etc.)
        cov = self.covmat_info.covmat
        if self.covmat_info.names==[spec.name for spec in self.spectra]:
            #Ordering is already ok
            return cov
        print "Covariance matrix is not in the same order as the 2pt measurement extensions...doing some damn fiddly"
        print "re-ordering, if I screw it up it's your fault for not putting your covariance in the right order"
        cov_starts = self.covmat_info.starts
        cov_lengths = self.covmat_info.lengths
        cov_names = self.covmat_info.names
        cov_ends=[cov_lengths[0]]
        for i in range(len(cov_names)-1):
            cov_ends.append(cov_ends[i]+cov_lengths[i+1])
        #print 'cov_lengths',cov_lengths
        #print 'cov_starts',cov_starts
        #print 'cov_ends',cov_ends
        assert cov_ends[-1]==cov.shape[0]

        total_l=0
        spec_inds=[]
        spec_names = [spec.name for spec in self.spectra]
        for spectrum in self.spectra:
            spec_inds.append(cov_names.index(spectrum.name))
            total_l+=cov_lengths[cov_names.index(spectrum.name)]
        cov_out=np.zeros((total_l,total_l))
        start_i=0
        #print spec_names
        #print spec_inds

        for ti,ind_i in zip(spec_names,spec_inds):
            start_j=0
            for tj,ind_j in zip(spec_names,spec_inds):
                cov_out[start_i:start_i+cov_lengths[ind_i],start_j:start_j+cov_lengths[ind_j]]=cov[cov_starts[ind_i]:cov_ends[ind_i],cov_starts[ind_j]:cov_ends[ind_j]]
                start_j+=cov_lengths[ind_j]
            start_i+=cov_lengths[ind_i]
        return cov_out

    def to_fits(self, filename, clobber=False):
        hdus = [fits.PrimaryHDU()]

        if self.covmat_info is not None:
            hdus.append(self.covmat_info.to_fits())

        for spectrum in self.spectra:
            if spectrum.windows not in window_types:
                raise NotImplementedError("Sorry - not yet coded general case with ell/theta window functions")
            hdus.append(spectrum.to_fits())

        if self.kernels is not None:
            for kernel in self.kernels:
                hdus.append(kernel.to_fits())

        hdulist = fits.HDUList(hdus)
        hdulist.writeto(filename, clobber=clobber)




    @classmethod
    def from_fits(cls, filename, covmat_name="COVMAT"):
        fitsfile = fits.open(filename)
        spectra = []
        kernels = {}
        windows = {}

        #Load the covariance matrix from the file, typically stored as
        #image data.
        #It turns out it is more conventient to load this first.
        #We can also use no covmat at all.
        if covmat_name is None:
            covmat_info = None
            covmat = None
        else:
            extension = fitsfile[covmat_name]
            covmat_info = CovarianceMatrixInfo.from_fits(extension)


        #First load all the spectra in the file
        #Spectra are indicated by the "2PTDATA" keyword
        #in the header being "T"
        for extension in fitsfile:
            if extension.header.get(TWOPOINT_SENTINEL):
                spectra.append(SpectrumMeasurement.from_fits(extension, covmat_info))

        #Each spectrum needs kernels, usually some n(z).
        #These were read from headers when we loaded the 2pt data above above.
        #Now we are loading those kernels into a dictionary
        for spectrum in spectra:
            for kernel in (spectrum.kernel1, spectrum.kernel2):
                if kernel not in kernels:
                    kernels[kernel] = NumberDensity.from_fits(fitsfile[kernel])

        #We might also require window functions W(ell) or W(theta). It's also possible
        #that we just use a single sample value of ell or theta instead.
        #If the spectra required it (according to the header) then we also
        #need to load those in.
        for spectrum in spectra:
            if spectrum.windows not in window_types and spectrum.windows not in windows:
                windows[spectrum.windows] = cls._windows_from_fits(fitsfile[windows])

        #return a new TwoPointFile object with the data we have just loaded
        return cls(spectra, kernels, windows, covmat_info)


                
    
    @classmethod
    def _windows_from_fits(cls, extension):
        raise NotImplementedError("non-sample window functions in ell/theta")
