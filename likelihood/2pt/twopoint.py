from astropy.io import fits
from astropy.table import Table
from enum34 import Enum
import numpy as np

#FITS header keyword indicating 2-point data
TWOPOINT_SENTINEL = "2PTDATA"
NZ_SENTINEL = "NZDATA"
COV_SENTINEL = "COVDATA"
PRECISION_SENTINEL = "PRECDATA"


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
    def __init__(self, name, bins, types, kernels, windows, angular_bin, value, angle=None, error=None):
        self.name = name
        self.bin1, self.bin2 = bins
        self.type1, self.type2 = types
        self.kernel1, self.kernel2 = kernels
        self.angular_bin = angular_bin
        self.angle = angle
        self.value = value
        self.windows = windows
        self.error = error

    def mask(self, mask):
        self.bin1 = self.bin1[mask]
        self.bin2 = self.bin2[mask]
        self.angular_bin = self.angular_bin[mask]
        self.angle = self.angle[mask]
        self.value = self.value[mask]

    def __len__(self):
        return len(self.value)

    def nbin(self):
        return np.max([self.bin1.max(), self.bin2.max()])

    def get_pair(self, bin1, bin2):
        w = (self.bin1==bin1) & (self.bin2==bin2)
        return self.angle[w], self.value[w]

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

        if windows != 'SAMPLE':
            raise NotImplementedError("Have not yet coded window functions for angular bins")

        #Now load the data
        data = extension.data
        bin1 = data['BIN1']
        bin2 = data['BIN2']
        angular_bin = data['ANGBIN']
        value = data['VALUE']
        angle = data['ANG'] if 'ANG' in data.names else None

        #Load a chunk of the covariance matrix too if present.
        if covmat_info is None:
            error = None
        else:
            error = covmat_info.get_error(name)

        return SpectrumMeasurement(name, (bin1, bin2), (type1, type2), (kernel1, kernel2), windows,
            angular_bin, value, angle, error)

    def to_fits(self):
        header = fits.Header()
        header[TWOPOINT_SENTINEL]= True
        header['EXTNAME']=self.name
        header['QUANT1'] = self.type1.value
        header['QUANT2'] = self.type2.value
        header['KERNEL_1'] = self.kernel1
        header['KERNEL_2'] = self.kernel2
        header['WINDOWS'] = 'SAMPLE' #NOT YET CODED ANYTHING ELSE

        columns = [
            fits.Column(name='BIN1', array=self.bin1, format='K'),
            fits.Column(name='BIN2', array=self.bin2, format='K'),
            fits.Column(name='ANGBIN', array=self.angular_bin, format='K'),
            fits.Column(name='VALUE', array=self.value, format='D'),
        ]
        if self.angle is not None:
            columns.append(fits.Column(name='ANG', array=self.angle, format='D'))

        extension = fits.BinTableHDU.from_columns(columns, header=header)
        return extension



class MatrixInfo(object):
    """Encapsulate a covariance or precision matrix and indices into it.

    Most of the behaviour of precision and covariance matrices is the same:
    we load, save, and index them in almost exactly the same way. So we define
    that behaviour in a base class and use two subclasses for the specifics.

    """
    sentinel = None
    def __init__(self, name, names, starts, lengths, matrix):
        super(MatrixInfo, self).__init__()
        self.name = name
        self.names = names
        self.starts = starts
        self.lengths = lengths
        self.matrix = matrix
        self.diagonal = matrix.diagonal()

    def to_fits(self):
        header = fits.Header()
        header[self.sentinel]= True
        header['EXTNAME'] = self.name
        for i,(start_index,name) in enumerate(zip(self.starts, self.names)):
            header['STRT_{}'.format(i)] = start_index
            header['NAME_{}'.format(i)] = name
        extension=fits.ImageHDU(data=self.matrix, header=header)
        return extension


    @classmethod
    def from_fits(cls, extension):
        if cls.sentinel is None:
            raise RuntimeError("""Do not use the MatrixInfo class directly; use a 
                subclass CovarianceMatrixInfo or PrecisionMatrixInfo.""")
        cov_name = extension.name
        covmat = extension.data
        header = extension.header
        assert cls.sentinel in header, """Tried to load a matrix of type {}
         from an extension {} but it does not contain that kind of data""".format(
            cls.sentinel, cov_name)
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
        for start, end in zip(start_indices[:-1], start_indices[1:]):
            lengths.append(end-start)
        if start_indices:
            lengths.append(covmat.shape[0]-lengths[-1])

        return cls(cov_name, measurement_names, start_indices, lengths, covmat)

class CovarianceMatrixInfo(MatrixInfo):
    sentinel=COV_SENTINEL
    def get_error(self, name):
        i = self.names.index(name)
        start = self.starts[i]
        end = start + self.lengths[i]
        return self.diagonal[start:end]**0.5

class PrecisionMatrixInfo(MatrixInfo):
    sentinel=PRECISION_SENTINEL



class TwoPointFile(object):
    def __init__(self, spectra, kernels, windows, covmat_info=None, precision_info=None):
        if windows is None:
            windows = {}
        self.spectra = spectra
        self.kernels = kernels
        self.windows = windows
        self.covmat_info = covmat_info
        self.precision_info = precision_info
        if covmat_info:
            self.covmat = covmat_info.matrix
        else:
            self.matrix = None
        if precision_info:
            self.precision = precision_info.matrix
        else:
            self.precision = None

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
            spectrum.mask(mask)
            print "Masking {} values in {}".format(mask.size-mask.sum(), spectrum.name)
            #record the mask vector as we will need it to mask the covmat
            masks.append(mask)
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
            spectrum.mask(mask)            
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

    def to_fits(self, filename, clobber=False):
        hdus = [fits.PrimaryHDU()]

        if self.covmat_info is not None:
            hdus.append(self.covmat_info.to_fits())

        if self.precision_info is not None:
            hdus.append(self.precision_info.to_fits())

        for spectrum in self.spectra:
            if spectrum.windows!="SAMPLE":
                raise NotImplementedError("Sorry - not yet coded general case with ell/theta window functions")
            hdus.append(spectrum.to_fits())

        for kernel in self.kernels.values():
            hdus.append(kernel.to_fits())

        hdulist = fits.HDUList(hdus)
        hdulist.writeto(filename, clobber=clobber)




    @classmethod
    def from_fits(cls, filename, covmat_name="COVMAT", precision_name=None):
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
        else:
            extension = fitsfile[covmat_name]
            covmat_info = CovarianceMatrixInfo.from_fits(extension)

        if precision_name is None:
            precision_info = None
            precision = None
        else:
            extension = fitsfile[precision_name]
            precision_info = PrecisionMatrixInfo.from_fits(extension)


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
            if spectrum.windows!="SAMPLE" and spectrum.windows not in windows:
                windows[spectrum.windows] = cls._windows_from_fits(fitsfile[windows])

        #return a new TwoPointFile object with the data we have just loaded
        return cls(spectra, kernels, windows, covmat_info, precision_info)


                
    
    @classmethod
    def _windows_from_fits(cls, extension):
        raise NotImplementedError("non-sample window functions in ell/theta")







