Franken Emu v1.0

Created by Earl Lawrence on 04/11/13

This code gives predictions for the nonlinear matter power spectrum
for a given set of cosmological parameters (flat wCDM). The k-range
and z-range are extended compared to the original Cosmic Emu and 
the hubble parameter can now be chosen freely if a slightly less
accurate answer is acceptable.

 The original emulation process is covered in "The Coyote Universe
III: Simulation Suite and Precision Emulator for the Nonlinear Matter
Power Spectrum" by Lawrence et al.  and the extension to high k and z
range is explained in "The Coyote Universe Extended: Precision
Emulation of the Matter Power Spectrum" by Heitmann et al.

In order to cover the larger ranges, a set of nested simulations
was added to the original Coyote Universe runs and power spectra
were created by combining different P(k)s from different simulation
boxes, leading us to call the new emulatore FrankenEmu.

The cosmological parameters and their ranges are:

0.120  < omega_m < 0.155 (including CDM and baryons)
0.0215 < omega_b < 0.0235
0.85   < n_s     < 1.05
0.61   < sigma_8 < 0.9
-1.30  < w       < -0.70
0      < z       < 1

The code will also ask at the start if the Hubble parameter should be
chosen internally based on CMB constraints or if the user wants
to provide a value. The range for the Hubble parameter is

55.0   < H0       < 85.0

The code will not produce an answer outside these ranges. The version 
that allows for the free choice of Hubble is slightly less accurate
as explained in detail in Heitmann et al.

The code works currently in C. A future release will include a Python
wrapper. The Gnu Scientific Library (GSL) is required for
compilation. In addition to the GSL, make sure you have all of

corrlengths.h
corrlengths_noh.h
design.h
design_noh.h
emu.c
emu_noh.c
hubble.c
kemu.h
ksim.h
meansd.h
meansd_noh.h
pcbasis.h
pcbasis_noh.h
pcweights.h
pcweights_noh.h
precisions.h
precisions_noh.h

A working program that demonstrates how to use everything is given
in main.c.   A makefile is provided that assumes the GSL is 
available and stored in a usual place.  If you have stored it 
somewhere else, then you're probably more qualified than any of us to fix
the compilation so that it works.  To compile the C example, type

$ make

and run "emu.exe" which will prompt you for an output file, parameter
values for each input, and the desired type of output.  The header
information contains the parameters and other information calculated
from the parameters.  The main output will be 2 columns. The first
column is k, the next column is the spectrum (in one of three forms)
for the chosen redshift z.

The main emulation function, emu, takes 3 parameters.  The first is a
vector with the cosmological parameters.  The second is a placeholder
vector which contains the k values and the spectra taken columnwise
from the format of the output file.  The third parameter is the output
type where 0 indicates the log(Delta^2 / k^(1.5)), 1 indicates
Delta^2, and 2 indicates P(k) = 2*(pi^2)*(Delta^2) / (k^3).

The function getH0fromCMB(double *x, double *stuff)
determines the H0 from the CMB constraint l_A=d_lss*pi/r_s=302.4 given
the cosmology specified in x.  This is the value used by the emulator.

Regarding the logo photo (CosmicEmu.jpg), the background comes from
NASA, ESA, and R. Massey at CalTech (
http://hubblesite.org/newscenter/archive/releases/2007/01/image/j/ ).
The emu is all over the Web and I don't know where it originated.  It
looks happy in space, I think.
