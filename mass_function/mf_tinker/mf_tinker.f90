!-----------------------------------------------------------------------------
!
! Halo multiplicity function in double precision, mf(lnnu), 
! where lnnu=ln(nu) and nu is called the "threshold," 
! given by nu = [delta_c/sigma(R,z)]^2. Here, delta_c=1.6865, sigma(R,z) 
! is the r.m.s. mass fluctuation within a top-hat smoothing of scale R 
! at a redshift z. 
!
! Ref. Tinker et al., ApJ, 688, 709 (2008)
! Note that "mf(lnnu)" here corresponds to f(sigma)/2 in this reference,
! where sigma=sqrt(1.6865/nu).
!
! IMPORTANT: Tinker et al.'s mass function is NOT normalized, i.e., 
!
! int_{-infinity}^{+infinity} dln(nu) mf(lnnu) = 1
!
! is NOT satisfied! Also, we use the fitting function for Delta=200.
! (That is, the mass is defined within a spherical region whose mean
! density is equal to 200 times the MEAN MASS density of the universe.)
!
! The halo mass function can be computed in the following way.
!
! (1) First, the comoving number density of halos per ln(nu) is related 
! to the mass function, dn/dM, as
!
! dM M dn/dM = rho_m dln(nu) mf(lnnu)
!
! where rho_m = (average mass density of the universe today)
!
! Dividing both sides by dM, one finds
!
! dn/dlnM = rho_m dln(nu)/dlnM mf(lnnu)/M
!
! (2) nu is more directly related to R, as nu=[delta_c/sigma(R,z)]^2.
! So, it is more convenient to write dn/dM as
!
! dn/dlnM = dlnR/dlnM dn/dlnR
!
! and compute dn/dlnR first, via
!
! dn/dlnR = rho_m dln(nu)/dlnR mf(lnnu)/M(R)
!
! Now, we can use M(R)=(4pi/3)rho_m R^3 to obtain
!
! dn/dlnR = (3/4pi) dln(nu)/dlnR mf(lnnu)/R^3
!
! with ln(nu) and dln(nu)/dlnR related to lnR via 
! - ln(nu) = 2*ln(1.6865)-ln[sigma^2(R,z)]
! - dln(nu)/dlnR = -dln[sigma^2(R)]/dlnR (which is indep. of z)
!
! (3) Once dn/dlnR is obtained as a function of R, one can 
! compute dn/dlnM using dlnR/dlnM=1/3, and obtain
!
! dn/dlnM = (1/3) dn/dlnR
!
! with M(R)=(4pi/3)rho_m R^3,
! and rho_m=2.775d11 (omega_matter h^2) M_solar Mpc^-3.
!
!
! << sample code >>
! 
! USE linearpk
! USE sigma
! double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
! double precision :: lnnu,dlnnudlnRh,Mh
! REAL :: chebev
! REAL :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
! character(len=128) :: filename
! integer :: n
! filename='wmap5baosn_max_likelihood_matterpower.dat'
! n=896 ! # of lines in the file
! CALL get_linearpk(filename,n)
! CALL compute_sigma2
! Rh=8. ! h^-1 Mpc
! lnRh = alog(Rh)
! if(lnRh>lnR2.or.lnRh<lnR1)then
!    print*,'lnR=',lnRh,' is out of range. It must be between',lnR1,' and',lnR2
!    stop
! else
!    lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
!    dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
!    lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
!    dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
!    dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0 ! in units of h^3 Mpc^-3
!    print*,'dn/dlnR at R=',Rh,' h^-1 Mpc is',dndlnRh,' h^3 Mpc^-3'
!    dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
!    Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
!    print*,'dn/dlnM at M=',Mh,' omega_matter h^-1 M_solar is (',dndlnMh,') h^3 Mpc^-3'
! endif
! end
!
! December 13, 2010
! E. Komatsu
!
!-----------------------------------------------------------------------------
MODULE mf_tinker
contains

DOUBLE PRECISION FUNCTION mf(lnnu,zz)
!!$ Halo multiplicity function per log(nu): 0.5*f(sigma), 
!!$ where sigma=deltac/sqrt(nu).
!!$ \int_{-infty}^infty mf(lnnu) dlnnu = 1 is NOT satisfied!
!!$ Ref: Eq.(3) of Tinker et al. (2008)
  IMPLICIT none
  DOUBLE PRECISION :: lnnu
  DOUBLE PRECISION :: nu,sigma, alpha, Delta=200.
  real(8):: zz
  double precision :: Ap,a,b
  double precision :: c=1.19d0 ! Delta=200
  Ap=0.186d0*(1.+ zz)**(-0.14)
  a=1.47d0*(1.+ zz)**(-0.06)
  alpha=10.00**(-1.0*(0.75/(dlog10(Delta/75.)))**1.2)
  b=2.57d0*(1.+ zz)**(-1.0*alpha)
  nu= exp(lnnu)
  sigma= 1.6865d0/dsqrt(nu)
  mf = 0.5d0*(Ap*((sigma/b)**(-a)+1d0)*dexp(-c/sigma**2d0))
  return
END FUNCTION mf
END MODULE mf_tinker
