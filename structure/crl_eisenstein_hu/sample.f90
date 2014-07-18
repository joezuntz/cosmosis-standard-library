MODULE sample
  IMPLICIT none

contains

  subroutine compute_pknowiggle(k,pk)



    k_ov_h=0.01d0 ! h Mpc^-1
    CALL eisensteinhu(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans)
    ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.002/Mpc
    ! Remember that ak is in units of h/Mpc whereas "k" in Eq.(74) is in units
    ! of 1/Mpc.
    pk=deltaR2*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0*growth_z0**2d0 &
         *trans**2d0*(k_ov_h*h0/0.002d0)**(ns-1d0) &
         *2d0*3.14159d0**2d0/k_ov_h**3d0
  end subroutine compute_pknowiggle

end MODULE sample
