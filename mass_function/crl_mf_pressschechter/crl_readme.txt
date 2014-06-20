*** Press&Schechter's Halo Mass Function ***
August 25, 2008: E.Komatsu

Ref. Press & Schechter, 187, 425 (1974) 
     See also Bond, Cole, Efstathiou & Kaiser, ApJ, 379, 440 (1991)
     See Sheth & Tormen, MNRAS, 329, 61 (2002) for the notation used here
     Note that "mf(lnnu)" here corresponds to nu*f(nu) in Sheth&Tormen(2002)

The halo mass function can be computed in the following way.

(1) First, the comoving number density of halos per ln(nu) is related 
to the mass function, dn/dM, as

  dM M dn/dM = rho_m dln(nu) mf(lnnu)

where 

- rho_m = (average mass density of the universe today)
- mf(lnnu) = (halo multiplicity function)
- lnnu=ln(nu) and nu is called the "threshold," 
given by nu = [delta_c/sigma(R,z)]^2, delta_c=1.6865, and sigma(R,z) 
is the r.m.s. mass fluctuation within a top-hat smoothing of scale R 
at a redshift z. 

Here, we provide a simple routine to compute the halo multiplicity
function, mf(lnnu), as a double-precision function. The argument, lnnu,
is also in double precision.

Dividing both sides by dM, one finds

  dn/dlnM = rho_m dln(nu)/dlnM mf(lnnu)/M

(2) nu is more directly related to R, as nu=[delta_c/sigma(R,z)]^2.
So, it is more convenient to write dn/dM as

  dn/dlnM = dlnR/dlnM dn/dlnR

and compute dn/dlnR first, via

  dn/dlnR = rho_m dln(nu)/dlnR mf(lnnu)/M(R)

Now, we can use M(R)=(4pi/3)rho_m R^3 to obtain

  dn/dlnR = (3/4pi) dln(nu)/dlnR mf(lnnu)/R^3

with ln(nu) and dln(nu)/dlnR related to lnR via 
- ln(nu) = 2*ln(1.6865)-ln[sigma^2(R,z)]
- dln(nu)/dlnR = -dln[sigma^2(R)]/dlnR (which is indep. of z)

(3) Once dn/dlnR is obtained as a function of R, one can 
compute dn/dlnM using dlnR/dlnM=1/3, and obtain

  dn/dlnM = (1/3) dn/dlnR

with M(R)=(4pi/3)rho_m R^3, 
and rho_m=2.775d11 (omega_matter h^2) M_solar Mpc^-3.

To compute sigma(R), it is necessary to use the input linear power spectrum.
We provide the sample data, "wmap5baosn_max_likelihood_matterpower.dat,"
which was generated using CAMB code for the maximum likelihood parameters
given in Table I of Komatsu et al.(2008) [WMAP 5-year interpretation paper] with 
"WMAP5+BAO+SN". The input file for CAMB is also provided 
(wmap5baosn_max_likelihood_params.ini). NOTE THAT THIS POWER SPECTRUM IS COMPUTED AT Z=0.

Another power spectrum, evolved back to z=30, is provided as 
"wmap5baosn_max_likelihood_matterpower_at_z=30.dat".

(1) Put "USE linearpk" at the beginning of the program 
(2) Put "USE sigma" at the beginning of the program
(3) CALL open_linearpk(filename,n) where "filename" contains the name of power spectrum 
file (character), and n is the number of lines in the file (integer)
(4) CALL compute_sigma2
(5) To compute ln(sigma^2), use a single-precision function CHEBEV(lnR1,lnR2,c,ndim,lnRh),
 where lnR1, lnR2, c and ndim have been pre-computed already by the step (3). 
The only argument you need to enter is the logarithm of the scale at which you want to 
compute ln(sigma^2), lnRh, which is in single precision. Note again that Rh is in units 
of h^-1 Mpc!
(6) To compute dln(sigma^2)/dlnRh, use a single-precision function 
CHEBEV(lnR1,lnR2,cder,ndim,lnRh), i.e., replace "c" in (5) with "cder".
(7) To convert ln(sigma^2) to ln(nu), use 2*ln(deltac)-ln(sigma^2), where deltac=1.6865
(8) To convert dln(sigma^2)/dlnRh to dln(nu)/dlnRh, use -dln(sigma^2)/dlnRh.
(9) Compute dn/dlnRh as (3/4pi)*dln(nu)/dlnRh*mf(lnnu)/Rh^3
(10) Compute dn/dlnMh as dndlnMh = dndlnRh/3 
(11) Convert Rh to Mh using Mh=(4pi/3)*2.775d11*Rh^3*omega_matter

- To compile and use the sample programs (given below), edit Makefile
and simply "./make"
- It will generate executables called "sample," "sample_mf_pressschechter" and
"compute_mf_pressschechter2"

=========================================================================
A simple program to use "mf_pressschechter.f90" (also included as "sample.f90" in
this directory):

 USE linearpk
 USE sigma
 double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
 double precision :: lnnu,dlnnudlnRh,Mh
 REAL :: chebev
 REAL :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
 character(len=128) :: filename
 integer :: n
 filename='wmap5baosn_max_likelihood_matterpower.dat'
 n=896 ! # of lines in the file
 CALL open_linearpk(filename,n)
 CALL compute_sigma2
 CALL close_linearpk
 Rh=8. ! h^-1 Mpc
 lnRh = alog(Rh)
 if(lnRh>lnR2.or.lnRh<lnR1)then
    print*,'lnR=',lnRh,' is out of range. It must be between',lnR1,' and',lnR2
    stop
 else
    lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
    dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
    lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
    dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
    dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0 ! in units of h^3 Mpc^-3
    print*,'dn/dlnR at R=',Rh,' h^-1 Mpc is',dndlnRh,' h^3 Mpc^-3'
    dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
    Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
    print*,'dn/dlnM at M=',Mh,' omega_matter h^-1 M_solar is (',dndlnMh,') h^3 Mpc^-3'
 endif
 end



Another program, "compute_mf_pressschechter.f90", for generating dn/dlnRh and dn/dlnMh 
for many R (h^-1 Mpc) and M (h^-1 M_solar) is

PROGRAM Compute_Mf_PressSchechter
  USE linearpk
  USE sigma
  ! A sample program for computing dn/dlnRh and dn/dlnMh using 
  ! Press&Schechter's formula.
  ! Units: 
  ! - Rh in h^-1 Mpc, 
  ! - Mh in omega_matter h^-1 M_solar, 
  ! - dndlnRh in h^3 Mpc^-3, and
  ! - dndlnMh in h^3 Mpc^-3
  ! August 25, 2008: E.Komatsu
  IMPLICIT none
  double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
  double precision :: lnnu,dlnnudlnRh,Mh
  real :: chebev
  real :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
  character(len=128) :: filename
  integer :: n
! read in and tabulate P(k)
  filename='wmap5baosn_max_likelihood_matterpower.dat'
  n=896 ! # of lines in the file
  CALL open_linearpk(filename,n)
! fit sigma^2(R) to Chebyshev polynomials
  CALL compute_sigma2
  CALL close_linearpk
! now output dn/dlnRh [in h^3 Mpc^-3] and dn/dlnMh [h^3 Mpc^-3] 
! as a function of R [h^-1 Mpc] and M [omega_matter h^-1 M_solar]...'
  open(1,file='Rh_dndlnRh.txt')
  open(2,file='Mhom0_dndlnMh.txt')
  do lnRh=lnR1,lnR2,0.01
     lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
     dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
     lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
     dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
     Rh=exp(lnRh) ! h^-1 Mpc
     dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0
     write(1,'(2E13.5)') Rh,dndlnRh
     dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
     Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
     write(2,'(2E13.5)') Mh,dndlnMh
  enddo
  close(1)
  close(2)
END PROGRAM Compute_Mf_PressSchechter

An extended version, "compute_mf_pressschechter2.f90", computes dn/dlnRh and dn/dlnMh 
at some redshift. The input power spectrum is given at z=30. 

PROGRAM Compute_Mf_PressSchechter2
  USE cosmo
  USE linearpk
  USE sigma
  USE growth
  ! A sample program for computing dn/dlnRh and dn/dlnMh using 
  ! Press&Schechter's formula at some redshift.
  ! Units: 
  ! - Rh in h^-1 Mpc, 
  ! - Mh in h^-1 M_solar, 
  ! - dndlnRh in h^3 Mpc^-3, and
  ! - dndlnMh in h^3 Mpc^-3
  ! August 25, 2008: E.Komatsu
  IMPLICIT none
  double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
  double precision :: lnnu,dlnnudlnRh,Mh
  real :: chebev
  real :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
  double precision :: zin,zout,g,scalefactor
  character(len=128) :: filename
  integer :: n
  external g
! Specify three cosmological parameters
! The data type has been defined in MODULE cosmo.
  ode0=0.723d0
  om0=0.277d0
  w=-1d0
! read in and tabulate P(k)
  filename='wmap5baosn_max_likelihood_matterpower_at_z=30.dat'
  n=896 ! # of lines in the file
  CALL open_linearpk(filename,n)
  zin=30d0
! fit sigma^2(R) to Chebyshev polynomials
  CALL compute_sigma2
  CALL close_linearpk
! ask for the output redshift
  print*,'output redshift?'
  read*,zout
! compute the growth factor
  CALL setup_growth
  scalefactor=g(zout)/g(zin)*(1d0+zin)/(1d0+zout)
  print*,'omega matter=',om0
  print*,'omega de=',ode0
  print*,'w=',w
! now output dn/dlnRh [in h^3 Mpc^-3] and dn/dlnMh [h^3 Mpc^-3] 
! as a function of R [h^-1 Mpc] and M [h^-1 M_solar] respectively...'
  open(1,file='Rh_dndlnRh.txt')
  open(2,file='Mh_dndlnMh.txt')
  do lnRh=lnR1,lnR2,0.01
     lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
     dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
     lnnu = 2d0*dlog(deltac)-dble(lnsigma2)-2d0*dlog(scalefactor) ! ln(nu)
     dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
     Rh=exp(lnRh) ! h^-1 Mpc
     dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0
     write(1,'(2E13.5)') Rh,dndlnRh
     dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
     Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
     write(2,'(2E13.5)') Mh*om0,dndlnMh
  enddo
  close(1)
  close(2)
END PROGRAM Compute_Mf_PressSchechter2
