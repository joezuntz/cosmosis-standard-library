!-----------------------------------------------------------------------------
!
! Compute sigma^2(R) from the linear matter power spectrum, and
! fit it by Chebyshev polynomials (single precision)
!
! << sample code >>
! 
! USE linearpk
! USE sigma
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
!    lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)
!    print*,'sigma(R) at R=',Rh,' h^-1 Mpc is',exp(0.5*lnsigma2)
!    dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh)
!    print*,'dln(sigma)/dlnR at R=',Rh,' h^-1 Mpc is',0.5*dlnsigma2dlnRh
! endif
! end
!
! August 23, 2008
! E. Komatsu
!
!-----------------------------------------------------------------------------
MODULE sigma
  INTEGER :: ndim=30
  REAL    :: c(30),cder(30)
  REAL    :: lnR1=-5.684 ! 0.0034Mpc/h, 1.8e4  solar mass
  REAL    :: lnR2=4.     ! 54.9Mpc/h, 7.5e16 solar mass
contains
  SUBROUTINE compute_sigma2
    REAL :: lnsigma2
    EXTERNAL lnsigma2
    CALL chebft(lnR1,lnR2,c,ndim,lnsigma2)
    CALL chder(lnR1,lnR2,c,cder,ndim)
    return
  END SUBROUTINE compute_sigma2
END MODULE sigma
!--------------------------------------------
REAL FUNCTION lnsigma2(lnR)
  IMPLICIT none
  REAL :: lnR
  DOUBLE PRECISION :: sigma2,r
  DOUBLE PRECISION :: ds2dk
  EXTERNAL ds2dk
  r = exp(lnR)
  CALL qromb(ds2dk,dlog(1d-4),dlog(2d1/r),sigma2,r)
  lnsigma2 = alog(REAL(sigma2))
  return
END FUNCTION lnsigma2
!--------------------------------------------
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
