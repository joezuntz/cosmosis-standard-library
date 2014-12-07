
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!subroutine Open_linearPk based on E. Komatsu's CRL routine which has been modified for use in CosmoSIS as we
!no longer read Pk from a file, all NR removed. EJ
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
MODULE LinearPk
  IMPLICIT none
  INTEGER :: jlo,ndata
  DOUBLE PRECISION, dimension(:), allocatable :: xa,ya,y2a
  DOUBLE PRECISION :: ak_min, ak_max
contains
  SUBROUTINE Open_LinearPk(aklin,pk_lin,n)
    IMPLICIT none
    integer, parameter :: dl=8
    real(dl),  dimension(:) :: aklin
    real(dl),  dimension(:) :: pk_lin
    integer n, ibcbeg,ibcend
    DOUBLE PRECISION :: ybcbeg,ybcend
    INTEGER :: i
    ndata=n
    ALLOCATE(xa(ndata),ya(ndata),y2a(ndata))
    xa=dlog(aklin)
    ya=dlog(pk_lin)
    !Natural spline settings
    ibcbeg = 2
    ybcbeg =0.0D+00
    ibcend = 2
    ybcend =0.0D+00
    CALL spline_cubic_set(ndata,xa,ya,ibcbeg,ybcbeg,ibcend,ybcend,y2a)
    ak_min=aklin(1)
    ak_max=aklin(ndata) 
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
  DOUBLE PRECISION :: a,b,h,x,y,ak,ypval, yppval
  if (ak .le. ak_min .or. ak .ge. ak_max) then
        linear_pk = 0.0
        return
  endif
  x  = dlog(ak)
  CALL spline_cubic_val( ndata, xa, ya, y2a, x, y, ypval, yppval )
  Linear_Pk = dexp(y)
  return
END FUNCTION Linear_Pk
