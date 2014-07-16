!Modified version of Komatsu's CRL function
!(http://www.mpa-garching.mpg.de/~komatsu/crl/)
!compute_pk_nowiggle for CosmoSIS - AM

! A sample program for computing the linear power spectrum
! of density fluctuations using Eisenstein&Hu's transfer function
! without the baryonic oscillation. See Eq.(29) of
! Eisenstein & Hu, ApJ, 496, 605 (1998)
! - k is in units of h Mpc^-1
! - P(k) is in units of h^-3 Mpc^3
! January 1, 2012: E.Komatsu


MODULE Compute_Pk_NoWiggle
  USE interface_tools
  USE pk_EisensteinHu
  IMPLICIT none

private
public compute_pknowiggle

contains


  subroutine compute_pknowiggle(k_ov_h,A_s,n_s,h0,om0,ob0,pk)
    !  k_ov_h in h Mpc^-1

    real(dl), INTENT(IN) :: k_ov_h,A_s,h0,om0,ob0,n_s
    real(dl),INTENT(INOUT) :: pk
    DOUBLE PRECISION :: trans
    CALL eisensteinhu(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans)
    ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.05/Mpc
    ! Remember that ak is in units of h/Mpc whereas "k" in Eq.(74) is in units
    ! of 1/Mpc. Be carefull you are using the same definition of As in other modules
    ! i.e same k_pivot

    pk=A_s*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0 &
         *trans**2d0*(k_ov_h*h0/0.05d0)**(n_s-1d0) &
         *2d0*3.14159265359d0**2d0/k_ov_h**3d0

    !  I commented growth_z0**2d0 so you have to put the z dependence by hands.
  end subroutine compute_pknowiggle


END MODULE Compute_Pk_NoWiggle
