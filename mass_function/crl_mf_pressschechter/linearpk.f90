!subroutine Open_linearPk has been slightly modified for use in CosmoSIS as we
!no longer read Pk from a file - EJ


!-----------------------------------------------------------------------------
!
! Read in, tabulate, and interpolate the linear matter power spectrum.
! The input file is in the following format (same as CAMB format):
!
! [1st column] k/h (i.e., wavenumber in units of h Mpc^-1) 
! [2nd column] h^3 P(k) (i.e., P(k) in units of h^-3 Mpc^3)
!
! << sample code >>
! 
! USE linearpk
! double precision :: linear_pk, k_ov_h
! character(len=128) :: filename
! integer :: n
! external linear_pk
! filename='wmap5baosn_max_likelihood_matterpower.dat'
! n=896 ! # of lines in the file
! CALL open_linearpk(filename,n)
! k_ov_h=1d0
! print*,'P(k) at k=',k_ov_h,' h Mpc^-1 is',linear_pk(k_ov_h),' h^3 Mpc^-3'
! end
!
! August 23, 2008
! E. Komatsu
!
!-----------------------------------------------------------------------------
MODULE LinearPk
  IMPLICIT none
  INTEGER :: jlo,ndata
  DOUBLE PRECISION, dimension(:), allocatable :: xa,ya,y2a
contains
 ! SUBROUTINE Open_LinearPk(filename,n)
  SUBROUTINE Open_LinearPk(aklin,pk_lin,n)
    IMPLICIT none
    integer, parameter :: dl=8
    real(dl),  dimension(:) :: aklin
    real(dl),  dimension(:) :: pk_lin
    integer n
    !DOUBLE PRECISION, dimension(:), allocatable :: aklin,pk_lin
    !CHARACTER(len=128) :: filename
    !INTEGER, intent(IN) :: n
    INTEGER :: i
    ndata=n
    ALLOCATE(xa(ndata),ya(ndata),y2a(ndata))
    !ALLOCATE(aklin(ndata),pk_lin(ndata))
    !open(1,file=filename,status='old',form='formatted')
    !print*,'read in '//trim(filename)
    !do i=1,ndata
     !  read(1,*)aklin(i),pk_lin(i) ! k/h & h^3 P(k) in CAMB format!
    !enddo
    !close(1)
    xa=dlog(aklin)
    ya=dlog(pk_lin)
    !DEALLOCATE(aklin,pk_lin)
    CALL spline(xa,ya,ndata,1.d30,1.d30,y2a)
    return
  END SUBROUTINE Open_LinearPk
  SUBROUTINE Close_LinearPk
    DEALLOCATE(xa,ya,y2a)
    return
  END SUBROUTINE Close_LinearPk
END MODULE LinearPk
DOUBLE PRECISION FUNCTION Linear_Pk(ak)
  Use LinearPk
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  x  = dlog(ak)
  CALL hunt(xa,ndata,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=a*ya(jlo)+b*ya(jlo+1)+((a**3-a)*y2a(jlo)+(b**3-b)*y2a(jlo+1))*(h**2)/6.
  Linear_Pk = dexp(y)
  return
END FUNCTION Linear_Pk
DOUBLE PRECISION FUNCTION dlnPkdlnk(ak)
  Use LinearPk
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  x  = dlog(ak)
  CALL hunt(xa,ndata,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=(ya(jlo+1)-ya(jlo))/h+(-(3.*a**2-1.)*y2a(jlo)+(3.*b**2-1.)*y2a(jlo+1))*h/6.
  dlnPkdlnk = y
  return
END FUNCTION DlnPkdlnk
