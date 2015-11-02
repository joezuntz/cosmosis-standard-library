# plc 2.0
(4 July 2015 version)

``plc`` is the public Planck Likelihood Code.  It provides C and Fortran 
libraries that allow users to compute the log likelihoods of the temperature, 
polarization, and lensing maps.  Optionally, it also provides a python version
of this library, as well as tools to modify the predetermined options for some 
likelihoods (e.g. changing the high-ell and low-ell lmin and lmax values of the
temperature).


## Installing

``plc`` can be installed either using the waf configuration and install tool, 
or using make.  Note that in the latter case, the installation procedure will
not test the availability of the tools and libraries needed by plc.


### Prerequisites

**Mandatory**

- C compiler, either gcc, clang or icc. 
- Fortran compiler, either ifort of gfortran
- ``blas/lapack`` library $
- ``cfitsio`` library $

**Optional**

- Python v>2.6, including the header and libraries.
- numpy $
- pyfits $
- cython $ 

All of the prerequisite labeled with a $ above can be installed by the waf
tool.

There are some incompatibilities between different versions of the C and
Fortran compilers.  In particular, gcc v4.9 and ifort v<14.1 are not
compatible.

Both ``blas/lapack`` and ``cfitsio`` need to be compiled as shared libraries
(option ``make --shared`` of the ``cfitsio`` install).

The optional prerequisites are only needed for the optional ``plc`` tools. 
They are not required for basic uses of the library.


### installing with waf

Waf (http://waf.io) is a tool based on python.  Installing with waf requires 
python v>2.6.  The tool first needs to be configured with
```
./waf configure [OPTIONS]
```
A complete list of the install options can be obtained by doing
```
./waf --help
```

The building and installation of the code is done using
```
./waf install
```


#### Simple configuration recipes

Detailed instructions on installing are found below, including a few simple and
typical cases.


##### MAC OS

Use the following instructions to install on Mac OS, letting the tool install
all of the optional prerequisites and using the stock blas/lapack library.  The
C compiler will be either gcc, icc or clang (tested in that order).  The
Fortran compiler will be either ifort or gfortran (tested in that order).
Installation will be performed in the current plc directory.
```
./waf configure --install_all_deps 
./waf install
```


##### Linux with mkl

Use the following instructions to install on a Linux machine with ifort and
mkl, letting the tool install all the optional prerequisites.  The mkl library
is available on the computer at the path ``$MKLROOT`` and its version is >=10.3 
(probably the case if your mkl library is post 2011).  The C compiler will be 
either gcc or clang (tested in that order). The C compiler will be either gcc, 
icc or clang (tested in that order).  The Fortran compiler will be eiher ifort
or gfortran (tested in that order).  Installation will be performed in the
current plc directory.
```
./waf configure --install_all_deps --lapack_mkl=$MKLROOT
./waf install
```


#### Detailed configuration instructions


##### Install the prerequisites automatically

```
./waf configure --install_all_deps [OPTIONS]
``` 
asks the tool to *try* to install all of the absent prerequisites.


##### Changing the install path

By default, plc is installed in the source directory.  To select a different
directory, use
```
./waf configure --prefix=/some/other/path [OPTIONS]
```


##### Selecting the C and Fortran compilers

By default, the tool will check for availability of gcc, icc, and clang, in
that order and select the first available compiler.  To force gcc, icc or clang
use one of ``--gcc``, ``--icc``, or ``--clang``.

By default, the tool will check for availability of ifort and gfortran, in that 
order and select the first available compiler.  To force ifort or gfortran, use 
one of ``--ifort`` or ``--gfortran``.


##### ``blas/lapack``

By default, on Mac OS, the stock ``blas/lapack`` will be used.  On Linux, the
waf tool  needs to be pointed to a particular library.  To use an alternative
library on Mac OS, or to point a particular library on Linux one can use the
following options.


###### ``mkl``

The simplest option is to use the Intel ``mkl`` library.  Assuming that the
library is installed at ``$MKLROOT`` and that the version of the library is
10.3 or higher (which should be the case if the library was installed post
2011), use
```
./waf configure --lapack_mkl=$MKLROOT [OPTIONS]
```
To use an earlier version of the library (10.0, 10.1, or 10.2), for example 
version 10.2 use
```
./waf configure --lapack_mkl=$MKLROOT --lapack_mkl_version=10.2 [OPTIONS]
```


##### Other ``blas/lapack``

One can provide the location of the blas/lapack install using either
1. ``--lapack_prefix=/some/path``, the include and library path will be 
``/some/path/lib`` and ``/some/path/include``, and the installer will try to 
link to ``libblas.so`` and ``liblapack.so``
2. ``--lapack_lib=/some/path/lib --lapack_include=/some/path/include``,
allowing the user to define the library and include path.  The installer will
try to link to ``libblas.so`` and ``liblapack.so``
3. ``--lapack_link='-I/some/path/include -L/some/path/lib -lmyblas -lmylapack',
allowing the user to define the full link line for the lapack/blas install.


##### Installing a basic blas/lapack automatically

Using either ``--install_all_deps`` or ``--install_lapack`` will ask the tool
to install a plain vanilla, non-optimized blas/lapack library.


#### ``cfitsio``

``plc`` requires a ``cftisio`` library compiled as a shared library (option 
``--shared`` when building cfitsio). The tool will look for the include and 
library in the usual locations.  If it fails, one can point it to a particular 
version of the ``cfitsio`` library using either
1. ``--cfitsio_prefix=/some/path``, the include and library will be 
``/some/path/include`` and ``/some/path/lib``
2. ``--cfitsio_lib=/some/path/lib --cfitsio_include=/some/path/include``,
allowing the user to define the library and include path. 

Alternatively, the ``cfitsio`` library can be installed by the tool using
either ``--install_all_deps`` or ``--install_cfitsio``.


### Installing with make

All prerequisites must be installed before running make.  The ``Makefile`` file 
must be adapted to the particular computer configuration.  The C and Fortran
compiler must be selected, and in the case of the ifort compiler also the list
of runtime libraries (needed to link with C).  Please look at the contents of
the file for details.

Once the ``Makefile`` file is correct, build and install with 
```
make install
```
The optional python tools can be built using
```
make python_install
```


### Testing the code

First source the ``clik_profile.[c]sh`` file to import the environment
variables, e.g. ``source source path/to/plc-2.0/bin/clik_profile.sh``.  Then,
running ``clik_example_C`` on any likelihood file should perform an automatic
test.  For example
```
$> clik_example_C plik_dx11dr2_HM_v18_TT.clik
----
clik version 14f84626be57 MAKEFILE
  smica
Checking likelihood 'plik_dx11dr2_HM_v18_TT.clik' on test data. got -380.979
expected -380.979 (diff -8.68062e-09)
----
...
```


## Basic use

The base installation of plc consists of a library (libclik.so) along with a C 
include and a Fortran module, which allows the user to
- initialize a likelihood
- inquire about the characteristics of a likelihood (list of spectra, lmin,
lmax and list of nuisance parameters)
- compute a log likelihood
- deallocate a likelihood


### Setting up the environment

To set the environment path and variables, source the ``clik_profile.sh`` (or 
``clik_profile.csh``) file.  This can be added to your ``.login`` or other 
shell initialization script.
The ``clik_profile.sh`` file is located in the ``bin`` directory of your plc 
install, i.e. in the newly created ``bin``  directory of the plc source 
directory if you have not selected a particular install path (see above).
Once the script is executed, the variable ``$CLIK_PATH`` will point to the root
directory of your plc install, and the plc ``bin`` directory will be included 
in your ``$PATH`` environment variable.


### Compiling with plc

The compilation and link options to compile a C code against the clik library
can be obtained using the script ``clik-config``.
Similarly, ``clik-config_f90`` provides the same information when compiling
a Fortran code.


### Clik file

The data for the library are stored in directories that will be named 
*clik files* in the following.  In these directories, one can find ascii
metadata files (``_mdb``) and data in multiple binary formats (FITS for
numerical data).


### Using the C Library

Include ``clik.h`` to use the plc library. 


#### Initialization

The function 
```
clik_object* clik_init(char* filepath, error **err);
``` 
initializes a CMB likelihood object of type ``clik_object`` from a clik file at 
location ``filepath``.  The error variable, ``err``, can be set to ``NULL``, in 
which case the function will exit in case of error.  Details on the optional 
error reporting can be found in ``errorlist.h``.

In most cases, more than one likelihood can be initialized in a code. 

Similarly,
```
clik_object* clik_lensing_init(char* filepath, error **err);
``` 
initializes a lensing likelihood.

During the initialization, an automatic test is performed and displayed on
stdout.


#### Obtaining the list of spectra, lrange, nuisance parameters


##### CMB likelihood

The function
```
void clik_get_lmax(clik_object *clikid, int lmax[6],error **err);
```
fills an array of 6 ``lmax`` integers.  Each element of the array corresponds to
the lmax of one of the CMB spectra needed by the CMB likelihood referred by
``clikid``.  The ordering is **TT EE BB TE TB EB**. The array elements are set
to -1 for each spectrum unused by the current likelihood.

As above, ``err`` can be set to ``NULL`` to ignore error management.  For
example, when used with the high-ell temperature likelihood, ``lmax`` will be 
set to `` 2508 -1 -1 -1 -1 -1``, meaning that only the TT spectrum is needed by 
the likelihood, and goes up to lmax = 2508.

The function 
```
int clik_get_extra_parameter_names(clik_object* clikid, parname **names, error
**err);
``` 
returns the number of nuisance parameters.  If the ``names`` argument is not
``NULL`` it will be filled with a pointer to a list of names of the nuisance
parameters.  After usage, the memory pointed to by ``names`` must be freed by
the user.  As above, ``err`` can be set to ``NULL`` to ignore error
management.


### Lensing likelihood

The function 
```
int clik_lensing_get_lmaxs(clik_lensing_object *lclik, int lmax[7],
error **err);
```
works like the CMB function ``clik_get_lmax``.
Note that the lmax array is longer than in the CMB case.  The order in this
case is **phiphi TT EE BB TE TB EB**.  In the future, this array may be
extended to include phi-CMB correlations. 

The function 
```
int clik_lensing_get_extra_parameter_names(clik_object* clikid, parname
**names, error **err);
``` 
works like the CMB one.


#### Computing a log likelihood

The function
```
double clik_compute(clik_object* clikid, double* cl_and_pars,error **err);
```
computes the log likelihood for the CMB likelihood pointed at by ``clikid`` and
for the Cls and nuisance parameters contained in the array ``cl_and_pars``.  As
above, ``err`` can be set to ``NULL`` to ignore error management.  The
``cl_and_pars`` array is an array of double precision values that contains the
CMB Cls and nuisance parameter.  The order in the array is as follows:

- first the Cls, from l=0 to lmax (inclusive) in the order **TT EE BB TE TB
EB**.
Note that this really does start from l=0, even although the effective lmin
of a given likelihood can be different from 0. Cls below the effective lmin 
will be ignored in the computation and can be set to any  value.
Only the spectra whose ``lmax`` (from the function ``clik_get_lmax``) is 
>-1 are present. 

- then the nuisance parameters in the order defined by 
``clik_get_extra_parameter_names``.

For example, for the low-ell TEB likelihood whose lmax = ``29 29 29 29 -1 -1`` 
and which has a single nuisance parameter, the Planck absolute calibration, the
``cl_and_pars`` array should be
```
c0_TT
c1_TT
c2_TT
...
c28_TT
c29_TT
c0_EE
c1_EE
...
c29_EE
c0_BB
...
c29_BB
c0_TE
...
c29_TE
planck_calibration
```

**Please note that this function expects Cl and not Dl=l(l+1)Cl/2pi.**
**Also note that the function returns a log likelihood and not a chi2.**

The function
```
double clik_lensing_compute(clik_object* clikid, double* cl_and_pars,error
**err);
```
provides the log likelihood for lensing, and works similarly to the CMB 
function.  The only difference is in the ordering of the ``cl_and_pars`` array, 
which must include the **phiphi** spectrum, before the CMB spectra.  Note that
CMB spectra are needed to compute the lensing normalization and biases.


#### Cleaning up

A CMB likelihood can be deallocated using 
```
void clik_cleanup(clik_object** pclikid);
``` 
and a lensing one with 
```
void clik_lensing_cleanup(clik_lensing_object **plclik);
```


#### Example code

The file in ``src/clik_example_C.c`` demonstrates the use of the C library.
This code compiles to an executable, which reads the path of a clik file as its
first argument and optionally Cl+nuisance ascii files as it next ones.
The likelihood is first initialized, and some information (lrange, nuisance
parameters) is displayed on screen.  Then each optional Cl+nuisance ascii file
is read and the log likelihood for each one is computed and displayed.

The Cl+nuisance files must be ascii files, containing one value on each line, 
forming a vector of values for the Cl (not Dl!) and nuisance parameters in the 
order expected for the compute functions. 


### Using the Fortran Library

Use the ``clik`` module in your code.


#### Initialization

The subroutine 
```
subroutine clik_init(clikid,filepath)
``` 
initializes a CMB likelihood with the data at ``filepath`` and fills ``clikid`` 
with an integer referring to the likelihood. 

Similarily 
```
subroutine clik_lensing_init(clikid,filepath)
``` 
initializes a lensing likelihood.

During the initialization, an automatic test is performed and displayed on
stdout.


#### Obtaining the list of spectra, lrange, nuisance parameters


##### CMB likelihood

The subroutine
```
subroutine clik_get_lmax(clikid,lmax)
```
fills an array of 6 ``lmax`` integers.  Each element of the array corresponds
to the lmax of one of the CMB spectra needed by the CMB likelihood referred to
by ``clikid``.  The ordering is **TT EE BB TE TB EB**.  The array elements are
set to -1 for each spectrum unused by the current likelihood.

For example, when used with the high-ell temperature likelihood, ``lmax`` will
be set to `` 2508 -1 -1 -1 -1 -1``, meaning that only the TT spectrum is needed
by the likelihood, and goes up to lmax = 2508.

The function 
```
integer(kind=4) function clik_get_extra_parameter_names(clikid,names)
``` 
returns the number of nuisance parameters. ``names`` is filled with a list of 
names of the nuisance parameters.

After usage, the memory pointed to by ``names`` must be freed by the user.


### Lensing likelihood

The subroutine
```
subroutine clik_lensing_get_lmax(clikid,lmax)
```
works like the CMB subroutine ``clik_get_lmax``.
Note that the lmax array is longer than in the CMB case.  The order in this
case is **phiphi TT EE BB TE TB EB**.  In the future, this array may be
extended to include phi-CMB correlations. 

The function 
```
integer(kind=4) function clik_lensing_get_extra_parameter_names(clikid,names)
``` 
works like the CMB one.


#### Computing a log likelihood

The function
```
real(kind=8) function clik_compute(clikid,cl_and_pars)
```
computes the log likelihood for the CMB likelihood referred to by ``clikid`` and
for the Cl and nuisance parameters contained in the array ``cl_and_pars``.
The ``cl_and_pars`` array is an array of doubles that contains the CMB Cls and 
nuisance parameters.  The order in the array is as follows:

- first the Cls, from l=0 to lmax (inclusive) in the order **TT EE BB TE TB
EB**. 
Note that this really does start from l=0, even although the effective lmin
of a given likelihood can be different from 0. Cls below the effective lmin 
will be ignored in the computation and can be set to any  value.
Only the spectra whose ``lmax`` (from the function ``clik_get_lmax``) is 
>-1 are present. 

- then the nuisance parameters in the order defined by 
``clik_get_extra_parameter_names``.

For example, for the low-ell TEB likelihood whose lmax = ``29 29 29 29 -1 -1`` 
and which has a unique nuisance parameter (the Planck absolute calibration),
the ``cl_and_pars`` array should be
```
c0_TT
c1_TT
c2_TT
...
c28_TT
c29_TT
c0_EE
c1_EE
...
c29_EE
c0_BB
...
c29_BB
c0_TE
...
c29_TE
planck_calibration
```

**Please note that the function expects Cl and not l(l+1)Cl/2pi.**
**Also note that the function returns a log likelihood and not a chi2.**

The function
```
real(kind=8) function clik_lensing_compute(clikid,cl_and_pars)
```
provides the log likelihood for a lensing map, and works similarly to the CMB 
function.  The only difference is in the ordering of the ``cl_and_pars`` array, 
which must include the **phiphi** spectrum, before the CMB spectra.  Note that
CMB spectra are needed to compute the lensing normalization and biases.


#### Cleaning up

A CMB likelihood can be deallocated using 
```
subroutine clik_cleanup(clikid)
``` 
and a lensing one with 
```
subroutine clik_lensing_cleanup(clikid)
```


#### Example code

The file in ``src/clik_example_F90.f90`` demonstrates the use of the C library.
This code compiles to an executable that reads the path of a clik file as its
first argument and optionally Cl+nuisance ascii files as it next ones.  The
likelihood is first initialized, and some information (lrange, nuisance
parameters) is displayed on screen.  Then each optional Cl+nuisance ascii file
is read and the log likelihood for each one is computed and displayed.

The Cl+nuisance files must be ascii files, containing one value on each line, 
forming a vector of values for the Cl (not Dl!) and nuisance parameters in the 
order expected for the compute functions. 


### Test tools

The executable 
```
clik_example_C /some/path/to/likelihood [Cl+nuisance.txt [Cl+nuisance.file
[...]]]
``` 
and 
```
clik_example_F90 /some/path/to/likelihood [Cl+nuisance.txt [Cl+nuisance.file
[...]]]
``` 
performs tests on a given likelihood file, displays some information and allows
the user to compute the log likelihood of Cl+nuisance vectors. 
The Cl+nuisance files must be ascii files, containing one value on each line, 
forming a vector of values for the Cl (not Dl!) and nuisance parameters in the 
order expected for the compute functions. 


## Optional tools

The optional tools are only available if the optional requirements are met.


### Using the python library

The library can be called from python by importing the ``clik`` python package.


#### Initialization

CMB Likelihoods are represented by an instance of the ``clik`` objects and
lensing likelihoods by an instance of ``clik_lensing`` objects.  Both are
initialized from a clik file.

```
import clik
CMBlkl = clik.clik("/some/path/to/clikfile")
lenslkl = clik.clik_lensing("/some/path/to/cliklensingfile")
```


#### Obtaining the list of spectra, lrange, nuisance parameters


##### CMB likelihood

The ``get_lmax`` method of the ``clik`` object returns a 6 element tuple 
containing the maximum multipole for each of the CMB spectra.  Please see the C
description above for further explanation.

The ``get_extra_parameter_names`` method of the ``clik`` object returns a tuple
containing the names of the nuisance parameters


##### Lensing likelihood

The ``get_lmax`` method of the ``clik_lensing`` object returns a 6 element
tuple containing the maximum multipole for each of the phi and CMB spectra.
Please see the C description above for further explanation.

The ``get_extra_parameter_names`` method of the ``clik_lensing`` object returns 
a tuple containing the names of the nuisance parameters

The ``get_clpp_fid`` method of the ``clik_lensing`` object returns the fiducial
lensing spectrum, used to perform the automatic test, and to compute the N0 and 
N1 biases.  Corrections to those biases are computed as perturbations around
this spectrum.  It returns an array that contains the phiphi lensing spectra
from 0 up to an lmax (inclusive) given by the first element of the tuple
returned by ``get_lmax``.

The ``get_clcmb_fid`` method of the ``clik_lensing`` object returns the
fiducial CMB spectra, used to perform the automatic test, and to compute the N0
and N1 biases.  Corrections to those biases are computed as perturbations
around those spectra.  It returns an array containing all of the CMB spectra
from 0 to the lmax given by the values at index >=1 in the tuple returned by
the ``get_lmax`` method.


#### Computing a log likelihood

The ``clik`` and ``clik_lensing`` objects are callable.  They expect a single 
argument ``cl_and_pars``, which must be either a 1-dimensional or a
2-dimensional array of floats.  The shape of the array must be either
``(ntot,)`` or ``(i,ntot)``.  If the array is 1-dimensional, it will be
reshaped internally to ``(1,ntot)``.

``ntot`` is the number of elements of the vector expected by the ``compute``
function described in the C and F90 versions of the library.  The function will
compute the log likelihood for each vector of spectra and nuisance parameters
and return an array of log likelihood values. 


#### Error handling

In the case where the library encounters an error, a ``clik.lkl.CError`` error
object will be raised.


### Example code

The script ``clik_example_py`` demonstrates how to use the python library.  Its
usage is similar to the ``clik_example_C`` and ``clik_exmaple_f90`` tools
described above.


### Displaying some imformation with ``clik_print``

The script ``clik_print`` displays on the screen some more information than the
``clik_example_XXX`` tools.  It can also be used to check the correct building
and installation of the library.

It expects a single argument on its command line, giving the path to a
likelihood file:
```
clik_print /where/is/my/clikfile
```


### Retrieving the selfcheck vector

The script `` clik_get_selfchek`` allows the user to retrieve the data used in
the automatic test performed at each initialization.  The scripts expects 2
arguments on its command line, the path to a clik file, and the name of a file
in which the Cl and parameters vector will be saved as an ascii file, in a
format suitable for the ``clik_example_XXX`` tools.  Indeed, calling any of the
``clik_example_XXX`` tools with the same likelihood and using the data vector
obtained this way should produce the same log likelihood value displayed during
the automatic test.
```
clik_get_selfcheck /some/clikfile /where/to/save/datavector
```


### Changing the lmin and lmax of some of the likelihoods

The ``clik_change_lrange`` allows users to change the lmin and lmax values of
the plik and commander likelihoods.  This allows for the reproduction of the
tests performed on the hybridization multipole, and on the high-ell lmax
described in the Planck 2015 Likelihood paper.
```
clik_change_lrange input_clik lmin lmax output_clik
```
The script expects an input clik file, the new lmin and lmax and an output clik 
file.  The input file will not be modified.  A value of -1 for lmin (or lmax)
directs the script not to modify the value.  The lmin and lmax are targets.
The output file will have lmin and lmax as close as possible to those values,
given the available data and the binning scheme (for plik).  They will be
displayed.  The script will fail if the path ``output_clik`` already exists.

Note that the output clik file will be stripped of the information needed by
the automatic test. When using this new clik file, no test will be performed
during initialization.

### Exploring the content of clik files

The scripts ``clik_cldf_ls`` and ``clik_cldf_dump`` allow users to explore the
contents of a clik file.  ``clik_cldf_ls`` allow users to browse the data tree
of the file, while ``clik_cldf_dump`` will display the content of a given
leaf. ``clik_cldf_ls`` works on the branches of the data tree (i.e. directories), 
while ``clik_cldf_dump`` only displays the content of the leafs (i.e. files).
Please note that when the leaf contains a long array, only a subset
of the first and last parts of the vector will be displayed.

The file can also be browsed in python using the ``clik.cldf`` module
```
import clik.cldf

# open file
clf = clik.cldf.File("plik_dx11dr2_HM_v18_TT.clik")

# get the content of the root
print clf.keys()

# is the path at clik a subtree or a leaf with data ?
# prints True in the former case 
print isinstance(clf["clik"],clik.cldf.File)

# get the list of content at path clik/lkl_0
print clf["clik/lkl_0"].keys()

# get the number of T channels in the likelihood
print  clf["clik/lkl_0/m_channel_T"]
```
