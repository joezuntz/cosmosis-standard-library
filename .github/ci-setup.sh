#The gnu science library
export GSL_INC=/usr/include
export GSL_LIB=/usr/lib

#The cfitsio FITS library
export CFITSIO_INC=/usr/include
export CFITSIO_LIB=/usr/lib

#The fftw3 Fourier transform library
export FFTW_LIBRARY=/usr/lib
export FFTW_INCLUDE_DIR=/usr/include

# BLAS and LAPACK are used from OpenBLAS, rather than
# from the macOS Accelerate framework. In this, we are
# following the lead of SciPy and NumPy.
export LAPACK_LIB=/usr/lib
export LAPACK_LINK="-L/usr/lib -llapack -lblas"
