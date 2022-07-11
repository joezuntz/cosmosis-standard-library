# EuclidEmulator2 (version 1.0.1)
This repository contains the source code of EuclidEmulator2, a fast and accurate tool to estimate the non-linear correction to the matter power spectrum.
In contrast to its predecessor EuclidEmulator, EuclidEmulator2 allows for 8-parameter cosmological models including massive neutrinos (assuming a degenerate hierarchy) and dynamical dark energy. EuclidEmulator2 is written in C++. For more information on EuclidEmulator please visit https://github.com/miknab/EuclidEmulator.

Authors:   M. Knabenhans, Pedro Carrilho<br/>
Date of last update:      March 2021<br/>
Reference: Euclid Consortium: Knabenhans et al. (submitted), <a>https://arxiv.org/abs/2010.11288</a><br/>

If you use EuclidEmulator2 in any way (for a publication or otherwise), please cite this paper.

<b>Contact information:</b> If you have any questions and/or remarks related to this work, please do not hesitate to send me an email (mischakATphysik.uzh.ch). For questions regarding the Python wrapper, contact Pedro Carrilho (p.gregoriocarrilhoATqmul.ac.uk).

## Currently implemented features
* emulation of the non-linear correction factor <i>B(k,z)</i>
* large allowed redshift interval: <i>z</i> in the interval [0.0,10.0]
* spatial scales spanning more than three orders of magnitude: 8.73 x 10<sup>-3</sup> <i>h</i> / Mpc ≤ <i>k</i> ≤ 9.41 <i>h</i> / Mpc.


* C++ executable
 + many different was to define the cosmologies:
   - direct definition of a single cosmology using command line parameters </li>
   - definition of several cosmologies through a parameter file </li>
   - definition of a cosmology through a CLASS or CAMB parameter file </li>
 + results are written to output file


* Python wrapper
  - Cosmology defined via parameter dictionary
  - Can output in custom k-range with extrapolation outside default range
  - Interfaces with class to compute full non-linear matter power spectrum

For a more extensive list of functionalities please the list of possible command line parameters shown [below](#building-and-installation). See the example jupyter notebook for details on running the python wrapper.

## Features yet to be implemented
* resolution correction emulator

## Quick start
### Prerequisites
In any case you need:
 * C++11 or later
 * GNU Scientific Library version 2.5 or higher (GSL; see https://www.gnu.org/software/gsl/)
 * g++ version 4.9.1 or higher

For the python wrapper, you need python with the following packages:
 * cython, numpy, scipy (essential)
 * jupyter, matplotlib (for plotting in the example notebook)
 * classy (for computing non-linear power spectrum only)

Note that the 3 essential packages will be installed during the wrapper installation via pip, if not present, but may not be recognised by anaconda later, which may lead to duplicate installations.

#### GSL install
On most machines, building GSL is relatively simple. To install it locally, e.g. in `~/local/gsl`, use
```
mkdir -p $HOME/local/gsl && cd $HOME/local
wget -c ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz -O - | tar zxv
```
The install procedure follows standard steps, but each one might take several minutes. Execute each command separately and only continue if there are no errors.
```
./configure --prefix=$HOME/local/gsl
make
make check
make install
```
 Once done, make sure to add the GSL library to your library path with
 ```
 export LD_LIBRARY_PATH=$HOME/local/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
 ```

### Test installations
The c++ code was successfully compiled on the following systems and environments:

* Mac OS X Mojave (10.14.6), with g++ version 9.3.0 and GSL version 2.6 (both g++ and GSL were installed with Homebrew)
* Linux (Red Hat 4.4.7-18), with g++ version 4.9.1 and GSL version 2.5
* Linux (Pop!_OS 20.04), with g++ version 9.3.0 and GSL version 2.6 (local GSL install)

### Get the code
If you have not done so already, either download this repository or clone it to your local host (under Linux you can get an unzipped tar-ball via
```
   wget -c https://github.com/PedroCarrilho/EuclidEmulator2/archive/pywrapper.tar.gz -O - | tar zxv
```

Alternatively, if you are only interested in the python wrapper, you can just install it using
```
pip install euclidemu2
```
See below for further installation instructions, in case pip install does not work.

### Building and installation
The current version of EuclidEmulator2 comes both as a command line interface (CLI) or a python wrapper.

#### CLI installation

In order to build the CLI, cd into the EuclidEmulator2 directory, check the `Makefile` and modify it as required (in particular the path of the `-I` and `-L` flag for GSL) and execute

```
   make
```
This will create an executable file called `ee2.exe`. In order to check if the build was successful, execute

```
   ./ee2.exe --help
```

You should now see the following output:

```
Highly accurate and efficient AI-based predictor of the non-linear correction of the matter power spectrum.
Usage:
  EuclidEmulator2 [OPTION...]

  -b, --Omega_b arg    baryon density parameter
  -m, --Omega_m arg    (total) matter density parameter
  -s, --Sum_m_nu arg   sum of neutrino masses in eV
  -n, --n_s arg        spectral index
  -H, --hubble arg     dimensionless Hubble parameter
  -W, --w_0 arg        time-independent dark energy equation-of-state
                       parameter
  -w, --w_a arg        time-dependent dark energy equation-of-state parameter
  -A, --A_s arg        spectral amplitude
  -z, --redshift arg   redshift at which the NLC shall be computed
  -i, --inifile arg    path to CLASS or CAMB parameter file
  -p, --parfile arg    path to EE2 parameter file
  -o, --outfile arg    output file name (default: nlc)
  -d, --outdir arg     output directory (default: results)
  -v, --verbosity arg  verbosity level (0, 1 or 2) (default: 0)
  -h, --help           Print usage
```

##### Sanity check
As a next step I recommend to go through the following:

  1. enter the EuclidEmulator2 directory (root level). Here you should see, amongst others, a `tests` directory.
  2. execute the command `mkdir results`
  3. execute the command `./tests/run_test.sh`

You should now see some output printed to screen that provides information about the cosmologies that are being treated by EuclidEmulator2. As a result, the following files should now be in the `results` directory:

```
  CAMB_example0.dat
  CLASS_example0.dat
  cmd_input_test.dat0.dat
  from_ee2_parfile0.dat
  from_ee2_parfile1.dat
  from_ee2_parfile2.dat
  from_ee2_parfile3.dat
```

In each file you should find a list of k modes (first column) and one or more B(k,z) column.


#### Python installation

To build the python wrapper, you can either use
```
pip install euclidemu2
```
to get the version from PyPI or clone this repo, cd into the EuclidEmulator2 directory and execute

```
pip install .
```

to install the version present on github.

##### Troubleshooting

Should these installation methods produce errors this can be due to two main issues:
* gsl path is not correct.
* g++ is not present or points to another compiler (such as clang on mac).

Solving the first issue is straightforward for a local installation by changing the paths to gsl present in setup.py, similarly to the CLI installation. The assumption for both ways of installing the python wrapper is that gsl is in `/usr/local`, such that the include files are in `/usr/local/include` and the lib files are in `/usr/local/lib`, so modify those lines in setup.py to the correct ones in your system and the installation should proceed.

Alternatively, you can pass additional arguments to pip to include the required paths. To install from PyPI, use instead
```
pip install --global-option=build_ext --global-option="-L/path/to/gsl/lib" --global-option="-I/path/to/gsl/include" euclidemu2
```
changing the `path/to/gsl` to your path.

To solve the second issue, install a new version of g++ and make sure that it is correctly linked, so that when you call `g++ --version` on a terminal you get an appropriate version of g++ and not clang.

To check if the installation was successful, just open python and do

```
import euclidemu2
```

A python notebook (test_euclid_emu2.ipynb) is also included, which includes an example for running the code and retrieving the boost factor. If it runs correctly, the installation was successful.

If you have further issues installing the python wrapper, don't hesitate to email Pedro Carrilho (p.gregoriocarrilhoATqmul.ac.uk) or post an issue on github.

### CLI Usage
There are several ways of calling EuclidEmulator2. The different options are briefly explained in the following:

* specify a cosmology directly by passing a value for each cosmological parameter (flags `-b, -m, -s, -n, -H, -W, -w, -A` and `-z` as well as there corresponding long versions). Notice that an arbitrary number of redshifts (`-z` flags) can be passed:

Example:
```
    ./ee2.exe -b 0.049 -m 0.319 -s 0.058 -n 0.96 -H 0.67 -W -1.0 -w 0.0 -A 2.1e-9\
              -z 0.0 -z 0.1081179 -z 0.4880429
```

* using the `-p` flag to pass an EuclidEmulator2 parameter file that contains at least one cosmology.

Example:
```
    ./ee2.exe -p tests/ee2_parfile.par
```

Checkout the file `tests/ee2_parfile.par` in order to learn more about how to structure such a parameter file.

* using the `-i` flag to pass a `CLASS` or `CAMB` parameter file.

Examples:
```
    ./ee2.exe -i tests/camb_parameter_file.ini -t CAMB
```
or
```
    ./ee2.exe -i tests/class_parameter_file.ini -t CLASS
```

Notice that the `-t` flag is mandatory in this case in order to tell the code whether it is reading a CLASS or a CAMB style \*.ini file.

### Python Usage

To run the python wrapper, first import the package via
```
import euclidemu2 as ee2
```
then create a python dictionary with the requested values of the cosmological parameters:
```
cosmo_par = {'As':2.1e-09, 'ns':0.966, 'Omb':0.04, 'Omm':0.3, 'h':0.68, 'mnu':0.15, 'w':-1.0, 'wa':0.0}
```
This can then be passed to the main function along with an array of redshifts
```
redshifts = [0,2,4,6,8,10]
k, b = ee2.get_boost(cosmo_par,redshifts)
```
resulting in an array `k` with the wavenumbers in units of <i>h</i> / Mpc and a dictionary `b` indexed with the same index as the redshift array, e.g. `b[1]` is an array corresponding to `redshifts[1]`, etc. This is always the case, even when only one redshift is requested, so accessing that array is always done via `b[0]`.

To calculate the boost for a custom range and sampling of k values, define an array with the requested values, e.g.
```
k_custom = np.geomspace(1e-4,100,1000)
```
and use
```
k, b = ee2.get_boost(cosmo_par,redshifts,k_custom)
```
where now `k=k_custom` and `b` gives the interpolated (extrapolated) boosts inside (outside) the default range.

If classy is installed, you can also get the full non-linear power spectrum via
```
k, pnl, plin, b = ee2.get_pnonlin(cosmo_par, redshifts, k_custom)
```
which will now output `pnl` in addition to the linear power spectrum `plin` and the boost `b`, which are all indexed in the same way as the boost from `get_boost`.

If classy is not installed, a warning will appear when loading `euclidemu2` and the `get_pnonlin` function will not work. The `get_boost` function will always work.

<b>Warning:</b> In the most recent versions of Python (e.g. 3.8) `classy` may not work unless it is the first package to be imported. This is taken into account when calling `euclidemu2`, but implies that `euclidemu2` must be the first package to be imported. This has been verified not to be a problem for older versions of python (e.g. 3.6).

See the python notebook (test_euclid_emu2.ipynb) for an example of a full run and more details.

## License
EuclidEmulator2 is free software, distributed under the GNU General Public License. This implies that you may freely distribute and copy the software. You may also modify it as you wish, and distribute these modified versions. You must always indicate prominently any changes you made in the original code and leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details.
