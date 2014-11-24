!Modified version of Komatsu's CRL function
!(http://www.mpa-garching.mpg.de/~komatsu/crl/) 
!compute_mf_pressschechter for CosmoSIS - EJ

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
                double precision :: lnsigma2,dlnsigma2dlnRh
                real :: Rh,lnRh
                real(dl) :: step
                integer :: int_lnRh
                integer :: n, nr,ii,nn
                real(dl),  dimension(:) :: k
                real(dl),  dimension(:) :: p
                type(massfunction) :: MassF
                double precision :: lnRh_1
                n=size(k)
                CALL Open_LinearPk(k,p,n)
                !fit sigma^2(R) to Chebyshev polynomials
                CALL compute_sigma2
                CALL close_linearpk

                step = 0.01
                nn= nint((lnR2 - lnR1)/step) !control variable for DO loop
                !allocate space for massfunction
                MassF%num_r = nr
                call allocate_mf(MassF)

                lnRh = lnR1
                DO ii=1,nn+1    
                !DO lnRh=lnR1,lnR2,0.1    
                        if(lnRh>lnR2) exit
                         lnRh_1 = (2.0D+00*lnRh - lnR1 - lnR2)/(lnR2 - lnR1)  !interval [-1,1]
                        CALL echebser0(DBLE(lnRh_1), c, ndim, lnsigma2)
                        CALL echebser0(DBLE(lnRh_1), cder, ndim, dlnsigma2dlnRh)
                        lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
                        dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
                        Rh=exp(lnRh) ! h^-1 Mpc
                        dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0
                        MassF%R_h(ii) = Rh
                        MassF%dn_dlnRh(ii) = dndlnRh
                        dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
                        Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
                        MassF%M_h(ii) = Mh
                        MassF%dn_dlnMh(ii) = dndlnMh
                        lnRh = lnR1 + ii*step
                END DO
        END SUBROUTINE compute_massfunction

END MODULE compute_mf_pressschechter
