# HMcode
HMcode readme
=============

This code is produces the matter power spectrum using the halo-model approach described in Mead et al. (2015). Appendix B of that paper details the methods for doing the calculation.

It should compile with any fortran compiler, and it doesn't need to be pointed to any libraries. I use 'ifort' and compile with '>ifort HMcode.f90'.

For the calculation a 'numin' and 'numax' value need to be set. These govern the range in nu that the 1-halo integral is taken over. The defaults I set are numin=0.3 and numax=5. In testing I looked at the power up to k=100. and numin=0.3 was sufficient to get convergence at this wavenumber for a standard cosmological model at all redshifts z<4. If one were interested in higher wavenumbers then numin would need to be decreased until convergence was achieved. numin should never be set to 0 because the mass function used in the code diverges at nu=0. numax=5 seems suitable for all practical purposes (note this means you are accounting for the effect of ~5 sigma fluctuations on the power-spectrum, which is quite generous).

When it starts the code fills up some arrays for k, z and power then calls 'assign_cosmology', which sets the cosmological parameters and tells the code where to look for an input linear power spectrum. The default is that this is taken from the Eistenstein + Hu (1998) approximation for the transfer function, but anything (e.g. a CAMB linear spectrum) could be wired in. See  the notes at the end of this README if you are interested in feeding HMcode a tabulated linear spectrum.

The code then loops through 'z' outer and 'k' producing power at each redshift and wave number. The ordering of loops (z then k) is because for each new redshift the halo-model calculation needs to call 'halomod_init' to fill up look-up tables, but then these are used for any 'k' at the redshift. At each new redshift the input linear spectrum is also renormalised (via the routine 'normalisation') so that it has the correct sigma8(z). After 'normalisation' a routine called 'fill_sigtab' is called which fills a look-up table of sigma(R) vs. R, which is useful in future calculations. This is a relatively expensive function to evaluate, hence the look-up table. 'halomod_init' then fills some look-up tables of various halo properties, such as mass, radius, nu, concentration etc. which are used in the one-halo integral.

Once these tables have been evaluated the halo-model integral can then be carried out. This calculation is done by the routine 'halomod', which calls 'p_1h' and 'p_2h' to evaluate 1- and 2-halo terms and then uses these to compute the full power spectrum. The power spectrum at each k and z is then added to an array which is printed out to power.dat (k, pow(z1), pow(z2), ...) when the code finishes. 

In testing I was able to get 16 redshifts, with 200 k-points, in 0.72 seconds (using ifort with -O3 optimisation). 

The file 'plot.p' is a gnuplot script to plot the output. It can be run using "gnuplot> load 'plot.p'".

Please let me know if you need any help running the code. Or if you have any questions whatsoever.

Alexander Mead
(am@roe.ac.uk)

######

Adding in a CAMB linear P(k)

Given the differences between CAMB and Eisenstein + Hu (1998) one might wish to make HMcode read in a linear CAMB power spectrum and work with that instead (or any other tabulated power spectrum). This is fine, and is what I did when comparing to Cosmic Emu in the paper (where the difference between CAMB and Eisenstein + Hu *was* important) but there is a subtle thing to bear in mind:

The halo-model calculation does integrals that are technically defined over the entire k range from 0 to inf. In practice contributions from very small and large scales are suppressed but the integration routines still need to know about the power at these scales sometimes, otherwise they may go bonkers. Obviously this is a problem given that one usually has a tabulated linear power spectrum defined on some finite range in k. The way I dealt with this was to read in the linear spectrum but then extrapolate if the code called 'p_lin' for k values below the minimum, or above the maximum, using physically motivated extrapolations. For example you know that the linear power spectrum is a power-law down to k->0 (\Delta^2 \propto k^{3+n}) and the high-k part can be approximated as \Delta^2 \propto k^{n-1}log(k)^2 . 

I have left my routines to do this in as 'find_Tk' and 'find_pk', and these carry out the correct extrapolation beyond the boundaries of either a P(k) table or T(k) table. These can be switched on using the 'itk' variable. Originally itk=3, which means the code uses Eisenstein + Hu (1998). If itk=4 is set then the code will look for an input table of k/h vs. T(k) and if itk=5 is set it will look for k/h vs. P(k). These input files need to be specified at run time (e.g. ./a.out input_tk.dat).


