MODULE pk_EisensteinHu
implicit none
private
public EisensteinHu
contains

SUBROUTINE EisensteinHu(ak,omegamh2,fb,T)
! REF: Eisenstein and Hu, ApJ, 496, 605 (1998), Eq.(29)-(31)
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: ak,omegamh2,fb ! 1/Mpc; omega_mh^2; omega_b/omega_m
  DOUBLE PRECISION, intent(OUT):: T
  DOUBLE PRECISION :: alpha,sound,shape,aq
  alpha= 1d0 &
       -0.328d0*dlog(431d0*omegamh2)*fb &
       +0.38d0*dlog(22.3d0*omegamh2)*fb**2d0
  sound= 44.5d0*dlog(9.83d0/omegamh2)  &
       /dsqrt(1d0+10d0*(fb*omegamh2)**(3d0/4d0))
  shape= omegamh2 &
       *(alpha+(1d0-alpha)/(1d0+(0.43d0*ak*sound)**4d0))
  aq= ak*(2.725d0/2.7d0)**2d0/shape
  T= dlog(2d0*dexp(1d0)+1.8d0*aq) &
       /(dlog(2d0*dexp(1d0)+1.8d0*aq) &
       +(14.2d0+731d0/(1d0+62.5d0*aq))*aq*aq)
  return
END SUBROUTINE EisensteinHu
end module pk_EisensteinHu