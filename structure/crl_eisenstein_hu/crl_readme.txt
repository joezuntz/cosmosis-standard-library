*** Linear matter density power spectrum from Eisenstein&Hu's transfer function without the baryonic oscillation ***
January 1, 2012: E.Komatsu

Here we provide a simple subroutine, eisensteinhu.f90, for computing Eisenstein&Hu's transfer function without the baryonic oscillation.

The linear matter power spectrum is given by [E.g., Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009)]

P(k,z) = Delta_R^2*(2*k^2/5/H0^2/Omega_M)^2
        *D(z)^2*T(k)^2*(k/k0)^(ns-1)
        *(2*pi^2)/k^3
where

- Delta_R^2 is the amplitude of the curvature perturbation in comoving gauge (or uniform curvature gauge) on superhorizon scale (e.g., Delta_R^2=2.46d-9)
- D(z) is the linear growth factor normalized such that (1+z)D(z)->1 during the matter era
- T(k) is the linear transfer function, computed using Eq.(29)-(31) of Eisenstein & Hu, ApJ, 496, 605 (1998)
- k0 is the wavenumber at which the initial power spectrum is 
normalized (e.g., k0=0.002/Mpc)

*** OUTPUT Of COMPUTE_PK_NOWIGGLE (wavenumber_pknowiggle.txt) ***
1st column: wavenumber [h Mpc^-1]
2nd column: power spectrum [h^-3 Mpc^3]

The output of this problem may be compared with the exact computation using the linear Boltzmann code. For this purpose we provide the sample data, "wmap5baosn_max_likelihood_matterpower.dat," which was generated using CAMB code for the maximum likelihood parameters given in Table I of Komatsu et al.(2008) [WMAP 5-year interpretation paper] with "WMAP5+BAO+SN". The input file for CAMB is also provided (wmap5baosn_max_likelihood_params.ini). NOTE THAT THIS POWER SPECTRUM IS COMPUTED AT Z=0.

Another power spectrum, evolved back to z=30, is provided as "wmap5baosn_max_likelihood_matterpower_at_z=30.dat".

- To compile and use the sample programs (given below), edit Makefile
and simply "make"
- It will generate executables called "sample" and "compute_pk_nowiggle"

=========================================================================
A simple program to use "eisensteinhu.f90" (also included as "sample.f90" in
this directory):

 double precision :: k_ov_h,pk,trans
 double precision :: om0=0.277d0,ob0=0.0459d0,h0=0.702d0
 double precision :: ns=0.962d0,deltaR2=2.46d-9
 double precision :: growth_z0=0.7646d0
 k_ov_h=0.01d0 ! h Mpc^-1
 CALL eisensteinhu(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans)
 ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.002/Mpc
 ! Remember that ak is in units of h/Mpc whereas "k" in Eq.(74) is in units 
 ! of 1/Mpc.
 pk=deltaR2*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0*growth_z0**2d0 &
   *trans**2d0*(k_ov_h*h0/0.002d0)**(ns-1d0) &
   *2d0*3.14159d0**2d0/k_ov_h**3d0
 print*,'P_11(k) at z=0 and k=',k_ov_h,' h Mpc^-1 is',pk,' h^-3 Mpc^3'
 end

A program for generating P(k) as a function of k, for a given redshift,
is provided as "compute_pk_nowiggle.f90". It looks like this:

PROGRAM Compute_Pk_NoWiggle
  USE cosmo
  USE growth
  ! A sample program for computing the linear power spectrum
  ! of density fluctuations using Eisenstein&Hu's transfer function
  ! without the baryonic oscillation. See Eq.(29) of 
  ! Eisenstein & Hu, ApJ, 496, 605 (1998)
  ! - k is in units of h Mpc^-1
  ! - P(k) is in units of h^-3 Mpc^3
  ! January 1, 2012: E.Komatsu
  IMPLICIT none
  integer :: j
  double precision :: g,z,D,trans
  double precision :: k_ov_h,pk,dlnk,lnk,zin
  double precision :: ob0,h0,ns,run,deltaR2
  external g
  character(len=128) :: filename
  integer :: n
! Specify three cosmological parameters
! The data type has been defined in MODULE cosmo.
  ode0=0.723d0
  om0=0.277d0
  w=-1d0
! Specify four more cosmological parameters
! These are not defined in MODULE cosmo.
  ob0=0.0459d0
  h0=0.702d0
  ns=0.962d0
  run=0d0
  deltaR2=2.46d-9
! tabulate g(z) by calling "setup_growth"
  call setup_growth
! ask for redshift
  print*,'redshift?'
  read*,z
  D=g(z)/(1d0+z) ! linear growth factor, normalized such that (1+z)D(z)=1 during the matter era
! now output P(k,z) as a function of k
  open(1,file='wavenumber_pknowiggle.txt')
  k_ov_h=1d-4 ! h/Mpc
  dlnk=2d-2
  do while (k_ov_h<6d3)
     CALL eisensteinhu(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans)
     ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.002/Mpc
     ! Remember that k_ov_h is in units of h/Mpc whereas "k" in Eq.(74) is in units 
     ! of 1/Mpc.
     pk=deltaR2*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0*D**2d0 &
       *trans**2d0*(k_ov_h*h0/0.002d0)**(ns-1d0+0.5d0*run*dlog(k_ov_h*h0/0.002d0)) &
       *2d0*3.14159d0**2d0/k_ov_h**3d0
     write(1,'(3E18.8)')k_ov_h,pk
     lnk=dlog(k_ov_h)+dlnk
     k_ov_h=dexp(lnk)
  enddo
  close(1)
END PROGRAM Compute_Pk_NoWiggle
