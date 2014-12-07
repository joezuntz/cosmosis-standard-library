!-----------------------------------------------------------------------------
!
! Compute sigma^2(R) from the linear matter power spectrum, and
! fit it by Chebyshev polynomials (single precision). EJ
!elise@fnal.gov
!-----------------------------------------------------------------------------
MODULE sigma
  INTEGER :: ndim=30
  DOUBLE PRECISION    :: c(30),cder(30)
  REAL    :: lnR1=-5.684 ! 0.0034Mpc/h, 1.8e4  solar mass
  REAL    :: lnR2=4.     ! 54.9Mpc/h, 7.5e16 solar mass
  DOUBLE PRECISION :: r
contains
  SUBROUTINE compute_sigma2
    DOUBLE PRECISION :: lnsigma2
    EXTERNAL lnsigma2
    CALL chebyshev_coefficients ( DBLE(lnR1), DBLE(lnR2), ndim, lnsigma2, c )
    CALL  dfrnt ( c, ndim, cder,DBLE(lnR1), DBLE(lnR2) ) !Modified Nov 2014 EJ
    return
  END SUBROUTINE compute_sigma2
END MODULE sigma
!--------------------------------------------
DOUBLE PRECISION FUNCTION lnsigma2(lnR)
  USE sigma
  IMPLICIT none
  DOUBLE PRECISION :: lnR
  REAL :: sigma2
  REAL :: wrap_ds2dk
  EXTERNAL wrap_ds2dk
  real abserr
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  integer, parameter :: key = 1
  integer n,ier
  integer neval
  r = DBLE(exp(lnR))
  call qag(wrap_ds2dk,REAL(dlog(1d-4)),REAL(dlog(2d1/r)), epsabs, epsrel, key, sigma2, abserr, neval, ier )
  lnsigma2 = DBLE(alog(sigma2))
  return
  END FUNCTION lnsigma2
!--------------------------------------------
!wrapper function as qag only works on REAL functions of one variable. EJ
REAL FUNCTION wrap_ds2dk(lnk)
  USE sigma
  DOUBLE PRECISION :: ds2dk
  REAL :: lnk
  EXTERNAL ds2dk
  wrap_ds2dk = REAL(ds2dk(DBLE(lnk),r))
  return
END FUNCTION wrap_ds2dk
DOUBLE PRECISION FUNCTION ds2dk(lnk,r)
  IMPLICIT none
  DOUBLE PRECISION :: lnk,k,r,x
  DOUBLE PRECISION :: linear_pk,win
  DOUBLE PRECISION :: pi= 3.14159265d0
  k = dexp(lnk)
  x = k*r
  win= (3d0/x)*(sin(x)/x**2d0-cos(x)/x)
  ds2dk = k**3d0*(linear_pk(k)*win**2d0)/(2d0*pi**2d0)
  return
END FUNCTION ds2dk
