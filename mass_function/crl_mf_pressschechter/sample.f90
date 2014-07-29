 USE linearpk
 USE sigma
 double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
 double precision :: lnnu,dlnnudlnRh,Mh
 REAL :: chebev
 REAL :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
 character(len=128) :: filename
 integer :: n
 filename='wmap5baosn_max_likelihood_matterpower.dat'
 n=896 ! # of lines in the file
 CALL open_linearpk(filename,n)
 CALL compute_sigma2
 CALL close_linearpk
 Rh=8. ! h^-1 Mpc
 lnRh = alog(Rh)
 if(lnRh>lnR2.or.lnRh<lnR1)then
    print*,'lnR=',lnRh,' is out of range. It must be between',lnR1,' and',lnR2
    stop
 else
    lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
    dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
    lnnu = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
    dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
    dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0 ! in units of h^3 Mpc^-3
    print*,'dn/dlnR at R=',Rh,' h^-1 Mpc is',dndlnRh,' h^3 Mpc^-3'
    dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
    Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
    print*,'dn/dlnM at M=',Mh,' omega_matter h^-1 M_solar is (',dndlnMh,') h^3 Mpc^-3'
 endif
 end
