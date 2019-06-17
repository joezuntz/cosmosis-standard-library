PROGRAM test
  IMPLICIT NONE
  REAL :: arr(5,2)
  REAL :: maxarr
  INTEGER :: sizearr(2)
  arr = reshape([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],[5,2])
  sizearr = shape(arr)
  maxarr = maxval(arr)
  WRITE (*,*) sizearr(1), sizearr(2)
  WRITE (*,*) maxarr
  CALL subrout(arr,sizearr)
END PROGRAM test

SUBROUTINE subrout(arr,sizearr)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sizearr(2)
  INTEGER :: i
  REAL, INTENT(IN) :: arr(sizearr(1),sizearr(2))
  DO i=1, sizearr(1)
     WRITE (*,*) arr(i,1), arr(i,2)
  END DO
END SUBROUTINE subrout
