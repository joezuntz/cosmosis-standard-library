PROGRAM Compute_Mf_ShethTormen2
  USE cosmo
  USE linearpk
  USE sigma
  USE growth
  ! A sample program for computing dn/dlnRh and dn/dlnMh using 
  ! Sheth&Tormen's formula at some redshift.
  ! Units: 
  ! - Rh in h^-1 Mpc, 
  ! - Mh in h^-1 M_solar, 
  ! - dndlnRh in h^3 Mpc^-3, and
  ! - dndlnMh in h^3 Mpc^-3
  ! August 25, 2008: E.Komatsu
  IMPLICIT none
  double precision :: deltac=1.6865d0,mf,dndlnRh,dndlnMh
  double precision :: lnnu,dlnnudlnRh,Mh
  real :: chebev
  real :: lnsigma2,lnRh,Rh,dlnsigma2dlnRh
  double precision :: zin,zout,g,scalefactor
  character(len=128) :: filename
  integer :: n
  external g
! Specify three cosmological parameters
! The data type has been defined in MODULE cosmo.
  ode0=0.723d0
  om0=0.277d0
  w=-1d0
! read in and tabulate P(k)
  filename='wmap5baosn_max_likelihood_matterpower_at_z=30.dat'
  n=896 ! # of lines in the file
  CALL open_linearpk(filename,n)
  zin=30d0
! fit sigma^2(R) to Chebyshev polynomials
  CALL compute_sigma2
  CALL close_linearpk
! ask for the output redshift
  print*,'output redshift?'
  read*,zout
! compute the growth factor
  CALL setup_growth
  scalefactor=g(zout)/g(zin)*(1d0+zin)/(1d0+zout)
  print*,'omega matter=',om0
  print*,'omega de=',ode0
  print*,'w=',w
! now output dn/dlnRh [in h^3 Mpc^-3] and dn/dlnMh [h^3 Mpc^-3] 
! as a function of R [h^-1 Mpc] and M [h^-1 M_solar] respectively...'
  open(1,file='Rh_dndlnRh.txt')
  open(2,file='Mh_dndlnMh.txt')
  do lnRh=lnR1,lnR2,0.01
     if(lnRh>lnR2) exit
     lnsigma2 = CHEBEV(lnR1,lnR2,c,ndim,lnRh)          ! ln(sigma^2)
     dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder,ndim,lnRh) ! dln(sigma^2)/dlnRh
     lnnu = 2d0*dlog(deltac)-dble(lnsigma2)-2d0*dlog(scalefactor) ! ln(nu)
     dlnnudlnRh = -dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
     Rh=exp(lnRh) ! h^-1 Mpc
     dndlnRh = (3d0/4d0/3.1415926535d0)*dlnnudlnRh*mf(lnnu)/dble(Rh)**3d0
     write(1,'(2E13.5)') Rh,dndlnRh
     dndlnMh = dndlnRh/3d0 ! in units of h^3 Mpc^-3
     Mh = (4d0*3.1415926535d0/3d0)*2.775d11*Rh**3d0 ! in units of omega_matter h^-1 M_solar
     write(2,'(2E13.5)') Mh*om0,dndlnMh
  enddo
  close(1)
  close(2)
END PROGRAM Compute_Mf_ShethTormen2
