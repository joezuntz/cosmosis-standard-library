!Modified version of Komatsu's CRL function
!(http://www.mpa-garching.mpg.de/~komatsu/crl/) 
!compute_mf_pressschechter for CosmoSIS - EJ


  ! A sample program for computing dn/dlnRh and dn/dlnMh using 
  ! Press&Schechter's formula.
  ! Units: 
  ! - Rh in h^-1 Mpc, 
  ! - Mh in omega_matter h^-1 M_solar, 
  ! - dndlnRh in h^3 Mpc^-3, and
  ! - dndlnMh in h^3 Mpc^-3
  ! August 25, 2008: E.Komatsu

MODULE compute_mf_pressschechter
  USE linearpk
  USE mf_pressschechter
  USE sigma
  USE interface_tools
  IMPLICIT none
  double precision :: deltac=1.6865d0,dndlnRh,dndlnMh

        contains

        SUBROUTINE compute_massfunction(k,p,MassF,nr)
                USE linearpk
                USE sigma
                USE interface_tools
                double precision :: lnnu,dlnnudlnRh,Mh
                real :: chebev
                real :: lnsigma2,Rh,lnRh,dlnsigma2dlnRh
                real(dl) :: step
                integer :: int_lnRh
                integer :: n, nr,ii,nn
                real(dl),  dimension(:) :: k
                real(dl),  dimension(:) :: p
                type(massfunction) :: MassF
                n=size(k)
                CALL Open_LinearPk(k,p,n)
                !fit sigma^2(R) to Chebyshev polynomials
                CALL compute_sigma2
                CALL close_linearpk

                ! now output dn/dlnRh [in h^3 Mpc^-3] and dn/dlnMh [h^3 Mpc^-3] 
                ! as a function of R [h^-1 Mpc] and M [omega_matter h^-1 M_solar]...'
               ! open(1,file='Rh_dndlnRh.txt')
                !open(2,file='Mhom0_dndlnMh.txt')
                step = 0.1
                nn= nint((lnR2 - lnR1)/step) !control variable for DO loop
                !allocate space for massfunction
                MassF%num_r = nr
                call allocate_mf(MassF)

                lnRh = lnR1
                DO ii=1,nn+1    
                !DO lnRh=lnR1,lnR2,0.1    
                        if(lnRh>lnR2) exit
                        lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
                        dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
                        lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
                        dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
                        Rh=exp(lnRh) ! h^-1 Mpc
                        dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0
                        MassF%R_h(ii) = Rh
                        MassF%dn_dlnRh(ii) = dndlnRh
                !        write(1,'(2E13.5)') Rh,dndlnRh
                        dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
                        Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
                        MassF%M_h(ii) = Mh
                        MassF%dn_dlnMh(ii) = dndlnMh
                !        write(2,'(2E13.5)') Mh,dndlnMh
                        lnRh = lnR1 + ii*step
                END DO
!  close(1)
!  close(2)
        END SUBROUTINE compute_massfunction

END MODULE compute_mf_pressschechter
