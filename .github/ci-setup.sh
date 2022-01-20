# This sets up various paths, most importantly $COSMOSIS_SRC_DIR
# which is used in most of the Makefiles to get the paths to cosmosis
# headers
source cosmosis-configure

export CC=gcc
export CXX=g++
export FC=gfortran

#The gnu science library
export GSL_INC=/usr/include
export GSL_LIB=/usr/lib/x86_64-linux-gnu

#The cfitsio FITS library
export CFITSIO_INC=/usr/include
export CFITSIO_LIB=/usr/lib/x86_64-linux-gnu

#The fftw3 Fourier transform library
export FFTW_INCLUDE_DIR=/usr/include
export FFTW_LIBRARY=/usr/lib/x86_64-linux-gnu

# BLAS and LAPACK are used from OpenBLAS, rather than
# from the macOS Accelerate framework. In this, we are
# following the lead of SciPy and NumPy.
export LAPACK_LIB=/usr/lib/
export LAPACK_LINK="-L/usr/lib/x86_64-linux-gnu -llapack -lblas"
