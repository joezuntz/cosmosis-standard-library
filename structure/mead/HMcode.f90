MODULE MHM


  TYPE cosmology
     REAL :: om_m, om_b, om_v, h, n, sig8, w, gamma
     REAL :: A
     REAL :: rnl, knl, sigv, neff
     INTEGER :: itk
     REAL, ALLOCATABLE :: rtab(:), sigtab(:)
     REAL, ALLOCATABLE :: ktab(:), tktab(:), pktab(:)

     !The two fitting parameters in the code
     real :: eta_0, As
  END TYPE cosmology

  TYPE tables
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
     INTEGER :: n
  END TYPE tables

  logical :: feedback !0 for no feedback, 1 for feedback
  REAL, PARAMETER :: pi=3.141592654

  CONTAINS

  subroutine meadfit_main()

  IMPLICIT NONE
  REAL :: p1h, p2h, pfull, plin
  REAL :: numin, numax, kmin, kmax, zmin, zmax
  CHARACTER(len=64) :: output
  REAL :: a, z
  INTEGER :: i, j, n, m, nk, nz
  TYPE(cosmology) :: cosi
  TYPE(tables) :: lut
  REAL, ALLOCATABLE :: k(:), simpower(:), ztab(:), ptab(:,:)
  character(64) :: input

  !Set to 1 to make halomodel calculation verbose
  feedback=.true.

  WRITE(*,*)
  WRITE(*,*) 'Welcome to Meadfit (must change name)'
  WRITE(*,*) '====================================='
  WRITE(*,*)

  !Set nu range for halomodel integration
  numin=0.1
  numax=5.

  !Set k range!
  nk=200
  kmin=0.001
  kmax=1.e4
  CALL fill_table(kmin,kmax,k,nk,1)

  WRITE(*,*) 'k min:', kmin
  WRITE(*,*) 'k max:', kmax
  WRITE(*,*) 'number of k:', nk
  WRITE(*,*)

  !Set the number of redshifts
  nz=16
  zmin=0.
  zmax=4.
  CALL fill_table(zmin,zmax,ztab,nz,0)

  WRITE(*,*) 'z min:', zmin
  WRITE(*,*) 'z max:', zmax
  WRITE(*,*) 'number of z:', nz
  WRITE(*,*)
  
  !Fill table for output power
  ALLOCATE(ptab(nz,nk))

  CALL assign_cosmology(cosi, input, 3)

  !Loop over redshifts
  DO j=1,nz

     !Sets the redshift
     z=ztab(j)

     !Initiliasation for the halomodel calcualtion
     !Also normalises power spectrum (via sigma_8)
     !and fills sigma(R) tables
     CALL halomod_init(z,numin,numax,lut,cosi)

     !Loop over k values
     DO i=1,SIZE(k)

        plin=p_lin(k(i),cosi)        

        CALL halomod(k(i),z,p1h,p2h,pfull,plin,lut,cosi)

        ptab(j,i)=pfull

     END DO

     IF(j==1) THEN
        WRITE(*,fmt='(A5,A7)') 'i', 'z'
        WRITE(*,fmt='(A13)') '   ============'
     END IF
     WRITE(*,fmt='(I5,F8.3)') j, ztab(j)

  END DO
  WRITE(*,*)

  output='power.dat'
  WRITE(*,fmt='(A19,A10)') 'Writing output to:', TRIM(output)
  WRITE(*,*)
  WRITE(*,*) 'The top row of the file contains the redshifts (the first entry is #####)'
  WRITE(*,*) 'Subsequent rows contain ''k'' and then the halo-model power for each redshift'
  OPEN(7,file=output)
  DO i=0,nk
     IF(i==0) THEN
        WRITE(7,fmt='(A20,40F20.10)') '#####', (ztab(j), j=1,nz)
     ELSE
        WRITE(7,fmt='(F20.10,40F20.10)') k(i), (ptab(j,i), j=1,nz)
     END IF
  END DO
  CLOSE(7)
  WRITE(*,*) 'Done'
  WRITE(*,*)

  end subroutine meadfit_main


  FUNCTION Delta_v(z,cosm)

    IMPLICIT NONE
    REAL :: Delta_v
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Virialised overdensity
    Delta_v=418.*(omega_m(z,cosm)**(-0.352))

  END FUNCTION Delta_v

  FUNCTION delta_c(z,cosm)

    IMPLICIT NONE
    REAL :: delta_c
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Linear collapse density
    delta_c=1.59+0.0314*log(sigma(8.,cosm))

  END FUNCTION delta_c

  FUNCTION eta(z,cosm)

    IMPLICIT NONE
    REAL :: eta
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !The first parameter here is 'eta_0' in Mead et al. (2015)
    !JAZ Moved this into the cosm type
    eta=cosm%eta_0-0.3*(sigma(8.,cosm))

  END FUNCTION eta

  FUNCTION kstar(z,cosm)

    IMPLICIT NONE
    REAL :: kstar
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !One-halo cut-off wavenumber
    kstar=0.584*(cosm%sigv)**(-1.)

  END FUNCTION kstar

  FUNCTION As(z,cosm)

    IMPLICIT NONE
    REAL :: As
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !This is the 'A' halo-concentration parameter in Mead et al. (2015)
    As=cosm%As

  END FUNCTION As

  FUNCTION fdamp(z,cosm)

    IMPLICIT NONE
    REAL ::fdamp
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Linear theory damping factor
    fdamp=0.188*(sigma(8.,cosm)**4.29)

  END FUNCTION fdamp

  FUNCTION alpha(z,cosm)

    IMPLICIT NONE
    REAL :: alpha
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !This uses the top-hat defined neff
    alpha=2.93*(1.77**cosm%neff)

  END FUNCTION alpha

  FUNCTION r_nl(lut)


    TYPE(tables), INTENT(IN) :: lut
    REAL :: r_nl

    !Calculates k_nl as 1/R where nu(R)=1.
    r_nl=exp(find(1.,lut%nu,log(lut%rr),3,3))

  END FUNCTION r_nl

  SUBROUTINE halomod(k,z,p1h,p2h,pfull,plin,lut,cosm)


    IMPLICIT NONE
    REAL, INTENT(OUT) :: p1h, p2h, pfull
    REAL, INTENT(IN) :: plin, k, z
    REAL :: c
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN
       p1h=0.
       p2h=0.
    ELSE
       p1h=p_1h(k,z,lut,cosm)
       p2h=p_2h(k,z,lut,plin,cosm)
    END IF

    c=alpha(z,cosm)
    pfull=(p2h**c+p1h**c)**(1./c)

  END SUBROUTINE halomod

  SUBROUTINE fill_table(min,max,arr,n,ilog)

    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL :: a, b
    REAL, ALLOCATABLE :: arr(:)
    INTEGER, INTENT(IN) :: ilog, n

    !Fills array 'arr' in equally spaced intervals
    !ilog=0 does linear spacing
    !ilog=1 does log spacing

    IF(ALLOCATED(arr)) DEALLOCATE(arr)

    ALLOCATE(arr(n))

    arr=0.

    IF(ilog==0) THEN
       a=min
       b=max
    ELSE IF(ilog==1) THEN
       a=log(min)
       b=log(max)
    END IF

    IF(n==1) THEN
       arr(1)=a
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=a+(b-a)*float(i-1)/float(n-1)
       END DO
    END IF

    IF(ilog==1) arr=exp(arr)

  END SUBROUTINE fill_table

  SUBROUTINE assign_cosmology(cosm, input, itk)


    IMPLICIT NONE
    CHARACTER(len=64) :: tk_file
    TYPE(cosmology) :: cosm
    LOGICAL :: lexist
    CHARACTER(len=64) :: input
    integer :: itk


    !itk==1 => Power-law models (i.e. T(k)=1.) (doesn't work with speed-ups)
    !itk==2 => DEFW transfer function
    !itk==3 => Eistenstein and Hu
    !itk==4 => Input CAMB T(k)
    !itk==5 => Input CAMB P(k)

    cosm%itk=itk

    cosm%om_m=0.3
    cosm%om_v=1.-cosm%om_m
    cosm%om_b=0.045
    cosm%h=0.7
    cosm%w=-1
    cosm%sig8=0.8
    cosm%n=0.97

    cosm%eta_0 = 0.603
    cosm%As = 3.13

    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_m:', cosm%om_m
    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_v:', cosm%om_v
    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_b:', cosm%om_b
    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'h:', cosm%h
    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'w:', cosm%w
    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'sig8:', cosm%sig8
    IF(feedback) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'n:', cosm%n   

    IF(cosm%itk==4) THEN

       !If the option 'itk' is switched to 4 then the code will look for an input k/h vs. T(k)
       !file. This must be speicifed when the program is run
       !(i.e. ./a.out input_tk.dat)

       CALL get_command_argument(1,input)
       IF(input=='') STOP 'ERROR: Please specify input T(k) file'
       INQUIRE(FILE=input, EXIST=lexist)  
       IF(lexist .EQV. .FALSE.) STOP 'ERROR: Specified T(k) file does not exist'

       CALL read_camb_tk(input,cosm)

    ELSE IF(cosm%itk==5) THEN

       !If the option 'itk' is switched to 4 then the code will look for an input k/h vs. P(k)
       !file. This must be speicifed when the program is run
       !(i.e. ./a.out input_pk.dat)

       CALL get_command_argument(1,input)
       IF(input=='') STOP 'ERROR: Please specify input P(k) file'
       INQUIRE(FILE=input, EXIST=lexist)  
       IF(lexist .EQV. .FALSE.) STOP 'ERROR: Specified P(k) file does not exist'

       CALL read_camb_pk(input,cosm)

    END IF

    IF(feedback) WRITE(*,*)    

  END SUBROUTINE assign_cosmology

  SUBROUTINE read_camb_tk(infile,cosm)


    IMPLICIT NONE
    INTEGER :: n, i
    CHARACTER(len=64) :: infile
    TYPE(cosmology) :: cosm

    !Reads-in a CAMB T(k) file

    n=file_length(infile)

    IF(ALLOCATED(cosm%ktab)) DEALLOCATE(cosm%ktab)
    IF(ALLOCATED(cosm%tktab)) DEALLOCATE(cosm%tktab)
    ALLOCATE(cosm%ktab(n),cosm%tktab(n))

    IF(feedback) WRITE(*,*) 'COSMOLOGY: Reading in CAMB transfer function' 
    IF(feedback) WRITE(*,*) 'COSMOLOGY: Length:', n
    IF(feedback) WRITE(*,*) 'COSMOLOGY: File name:', infile
    OPEN(7,file=infile)
    DO i=1,n
       READ(7,*) cosm%ktab(i), cosm%tktab(i)
    END DO
    CLOSE(7)
    IF(feedback) WRITE(*,*) 'COSMOLOGY: Finished'

    !Normalise so that T(k<<k_eq)=1
    cosm%tktab=cosm%tktab/cosm%tktab(1)

  END SUBROUTINE read_camb_tk

  SUBROUTINE read_camb_pk(infile,cosm)


    IMPLICIT NONE
    INTEGER :: n, i
    REAL :: c
    REAL, ALLOCATABLE :: ksmall(:), dksmall(:)
    CHARACTER(len=64) :: infile
    TYPE(cosmology) :: cosm
    REAL, PARAMETER :: pi=3.141592654

    !Reads-in a CAMB P(k) file

    n=file_length(infile)

    IF(ALLOCATED(cosm%ktab)) DEALLOCATE(cosm%ktab)
    IF(ALLOCATED(cosm%pktab)) DEALLOCATE(cosm%pktab)
    ALLOCATE(cosm%ktab(n),cosm%pktab(n))

    IF(feedback) WRITE(*,*) 'COSMOLOGY: Reading in CAMB power spectrum, length:', n
    IF(feedback) WRITE(*,*) 'COSMOLOGY: File name:', infile
    OPEN(32,file=infile)
    DO i=1,n
       READ(32,*) cosm%ktab(i), cosm%pktab(i)
    END DO
    CLOSE(32)
    IF(feedback) WRITE(*,*) 'COSMOLOGY: Finished'

    !Convert P(k) to \Delta^2(k)
    cosm%pktab=cosm%pktab*(cosm%ktab**3.)/(2.*(pi**2.))

  END SUBROUTINE read_camb_pk

  SUBROUTINE normalisation(z,cosm)


    IMPLICIT NONE
    REAL :: sigi, growz
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm

    !This normalises the power spectrum via sigma8

    !Must be set for the rest of the calculation to work properly
    cosm%A=1.

    !These tables *must* be deallocated for the normalisation calculation to work properly
    !This means they should be re-allocated afterwards (using the fill_sigtab routine)
    IF(ALLOCATED(cosm%rtab)) DEALLOCATE(cosm%rtab)   
    IF(ALLOCATED(cosm%sigtab)) DEALLOCATE(cosm%sigtab)
    
    sigi=sigma(8.,cosm)

    IF(feedback) THEN
       WRITE(*,*) 'NORMALISATION: Initial sigma8:', sigi
    END IF

    growz=grow(z,cosm)
    cosm%A=growz*cosm%sig8/sigi

    sigi=sigma(8.,cosm)

    IF(feedback) THEN
       WRITE(*,*) 'NORMALISATION: Growth factor:', growz
       WRITE(*,*) 'NORMALISATION: Normalisation factor:', cosm%A
       WRITE(*,*) 'NORMALISATION: Target sigma8 at z=0:', cosm%sig8
       WRITE(*,*) 'NORMALISATION: Final sigma8 (calculated at z):', sigi
       WRITE(*,*) 'NORMALISATION: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE normalisation

  SUBROUTINE allocate_LUT(lut)


    TYPE(tables) :: lut
    INTEGER :: n

    !Allocates memory for the look-up tables

    n=lut%n

    ALLOCATE(lut%zc(n),lut%m(n),lut%c(n),lut%rv(n))
    ALLOCATE(lut%nu(n),lut%rr(n),lut%sigf(n),lut%sig(n))

    lut%zc=0.
    lut%m=0.
    lut%c=0.
    lut%rv=0.
    lut%nu=0.
    lut%rr=0.
    lut%sigf=0.
    lut%sig=0.

  END SUBROUTINE allocate_LUT

  SUBROUTINE deallocate_LUT(lut)


    TYPE(tables) :: lut

    !Deallocates look-up tables

    DEALLOCATE(lut%zc,lut%m,lut%c,lut%rv,lut%nu,lut%rr,lut%sigf,lut%sig)

  END SUBROUTINE deallocate_LUT

  SUBROUTINE halomod_init(z,nu_min,nu_max,lut,cosm)


    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    INTEGER :: i, imin, imax, n
    REAL :: rrmin, rrmax, dr, nu_min, nu_max, Dv, dc, ms, f
    REAL*8 :: rknl, rneff, rncur
    REAL, ALLOCATABLE :: rg_up(:), rg_dn(:), nu_up(:), nu_dn(:), mg_up(:), mg_dn(:)
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut

    !Halo-model initialisation routine
    !First normalises the linear spectrum to have the correct sigma8
    !Then fills sigma(R) look-up tables
    !The computes other tables necessary for the one-halo integral

    !This then normalises the power spectrum at your z
    CALL normalisation(z,cosm)

    !Now fill the sigma tables
    CALL fill_sigtab(cosm)

    !Find value of sigma_v
    cosm%sigv=sqrt(dispint(cosm))

    n=1000

    ALLOCATE(rg_up(n),rg_dn(n),nu_up(n),nu_dn(n),mg_up(n),mg_dn(n))

    IF(feedback) WRITE(*,*) 'HALOMOD: Filling look-up tables'
    IF(feedback) WRITE(*,*) 'HALOMOD: Tables being filled at redshift:', z       
    IF(feedback) WRITE(*,*) 'HALOMOD: sigv [Mpc/h]:', cosm%sigv

    IF(ALLOCATED(lut%rr)) CALL deallocate_LUT(lut)

    dc=delta_c(z,cosm)

    !It seems this number can be quite large without changing P(k) much.
    !This governs the accuracy (sampling of log r) dr=0.05 seems okay
    !dr=0.02 seems more conservative and does not change run-time so much
    dr=0.02

    DO i=1,n

       !rg_up will start at r=1 Mpc/h
       rg_up(i)=10.**(dr*float(i-1))
       mg_up(i)=mass_r(rg_up(i),cosm)
       nu_up(i)=dc/sigma(rg_up(i),cosm)

       !Finishes once the numax value has been reached
       IF(nu_up(i)>nu_max) THEN
          imax=i
          EXIT
       END IF

    END DO

    DO i=1,n

       rg_dn(i)=10.**(-dr*float(i))
       mg_dn(i)=mass_r(rg_dn(i),cosm)
       nu_dn(i)=dc/sigma(rg_dn(i),cosm)

       !Finishes if numin is reached *or* of the halo mass is less than *1/h* Solar mass (!!)
       IF(nu_dn(i)<nu_min .OR. mg_dn(i)<1.) THEN
          imin=i
          EXIT
       END IF

    END DO

    lut%n=imin+imax

    IF(feedback) WRITE(*,*) 'HALOMOD: Entries in look up table:', lut%n

    CALL allocate_LUT(lut)

    DO i=1,imin
       lut%rr(i)=rg_dn(imin-(i-1))
       lut%nu(i)=nu_dn(imin-(i-1))
       lut%m(i)=mg_dn(imin-(i-1))
    END DO

    DO i=1,imax
       lut%rr(i+imin)=rg_up(i)
       lut%nu(i+imin)=nu_up(i)
       lut%m(i+imin)=mg_up(i)
    END DO

    lut%sig=dc/lut%nu

    IF(feedback) WRITE(*,*) 'HALOMOD: m, r, nu tables filled'

    !Fills up a table for sigma(fM) for Bullock c(m) relation
    !This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)
    f=0.01**(1./3.)
    DO i=1,lut%n
       lut%sigf(i)=sigma(lut%rr(i)*f,cosm)
    END DO
    IF(feedback) WRITE(*,*) 'HALOMOD: sigf tables filled'  

    !Fill virial radius table using real radius table
    Dv=Delta_v(z,cosm)
    lut%rv=lut%rr/(Dv**(1./3.))

    IF(feedback) WRITE(*,*) 'HALOMOD: rv tables filled'  
    IF(feedback) WRITE(*,*) 'HALOMOD: nu min:', lut%nu(1), 'target:', nu_min
    IF(feedback) WRITE(*,*) 'HALOMOD: nu max:', lut%nu(lut%n),'target:',  nu_max
    IF(feedback) WRITE(*,*) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
    IF(feedback) WRITE(*,*) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)
    IF(feedback) WRITE(*,*) 'HALOMOD: M min [Msun/h]:', lut%m(1)
    IF(feedback) WRITE(*,*) 'HALOMOD: M max [Msun/h]:', lut%m(lut%n)

    !Find non-linear radius and scale
    cosm%rnl=r_nl(lut)
    cosm%knl=1./cosm%rnl

    !Numerical differentiation to find effective index at collapse
    cosm%neff=-3.-derivative_table(log(cosm%rnl),log(lut%rr),log(lut%sig**2.),3,3)

    IF(feedback) WRITE(*,*) 'HALOMOD: r_nl [Mpc/h]:', cosm%rnl
    IF(feedback) WRITE(*,*) 'HALOMOD: k_nl [h/Mpc]:', cosm%knl
    IF(feedback) WRITE(*,*) 'HALOMOD: n_eff:', cosm%neff

    CALL conc_bull(z,cosm,lut)

    IF(feedback) WRITE(*,*) 'HALOMOD: c tables filled'
    IF(feedback) WRITE(*,*) 'HALOMOD: c min [Msun/h]:', lut%c(lut%n)
    IF(feedback) WRITE(*,*) 'HALOMOD: c max [Msun/h]:', lut%c(1)
    IF(feedback) WRITE(*,*) 'HALOMOD: Done'

    IF(feedback) WRITE(*,*)

    feedback=.false.

  END SUBROUTINE halomod_init

  SUBROUTINE conc_bull(z,cosm,lut)


    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut
    REAL :: A, ch, be, zf, g_lcdm, g_wcdm, w
    integer i
    !Calculates the Bullock et al. (2001) mass-concentration relation

    A=As(z,cosm)

    !Fill the collapse z look-up table
    CALL zcoll_bull(z,cosm,lut)

    !Fill the concentration look-up table
    DO i=1,lut%n

       zf=lut%zc(i)
       lut%c(i)=A*(1.+zf)/(1.+z)

       !Dolag2004 prescription for adding DE dependence
       IF(cosm%w .NE. -1.) THEN

          zf=100.

          g_wcdm=grow(zf,cosm)

          w=cosm%w
          cosm%w=-1.
          g_lcdm=grow(zf,cosm)
          cosm%w=w

          lut%c(i)=lut%c(i)*(g_wcdm/g_lcdm)

       END IF

    END DO

  END SUBROUTINE conc_bull

  SUBROUTINE zcoll_bull(z,cosm,lut)


    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut
    REAL :: dc
    REAL :: amin, amax, af, zf, RHS
    REAL, ALLOCATABLE :: af_tab(:), grow_tab(:)
    INTEGER :: i, j, ntab

    !This fills up the halo formation redshift table as per Bullock relations

    !Needs to interpolate g(z) which should be pretty linear for a<0.05 
    !in 'g(a) vs a' space for all standard cosmologies
    amin=0.05
    amax=1./(1.+z)

    ntab=64

    ALLOCATE(af_tab(ntab),grow_tab(ntab))

    af_tab=0.
    grow_tab=0.

    !Fill the growth function tables once only
    DO j=1,ntab

       af_tab(j)=amin+(amax-amin)*float(j-1)/float(ntab-1)
       zf=-1.+1./af_tab(j)
       grow_tab(j)=grow(zf,cosm)

    END DO

    dc=delta_c(z,cosm)

    !Do numerical inversion
    DO i=1,lut%n

       RHS=dc*grow(z,cosm)/lut%sigf(i)

       IF(RHS>grow_tab(ntab)) THEN
          !This is the case of 'halo forms in the future'
          !in this case set formation redshift to current redshift
          zf=z
       ELSE
          af=find(RHS,grow_tab,af_tab,3,3)
          zf=-1.+1./af
       END IF

       lut%zc(i)=zf

    END DO

    DEALLOCATE(af_tab,grow_tab)

  END SUBROUTINE zcoll_bull

  FUNCTION mass_r(r,cosm)

    IMPLICIT NONE
    REAL :: mass_r, r
    TYPE(cosmology) :: cosm
    REAL, PARAMETER :: pi=3.141592654

    !Relation between mean cosmological mass and radius

    mass_r=(4.*pi/3.)*cosmic_density(cosm)*(r**3.)

  END FUNCTION mass_r

  FUNCTION cosmic_density(cosm)

    IMPLICIT NONE
    REAL :: cosmic_density
    TYPE(cosmology) :: cosm

    !Cosmological density at z=0 in M_sun per Mpc^3 with h factors included. 
    !The constant does this.

    cosmic_density=(2.775e11)*cosm%om_m

  END FUNCTION cosmic_density

  FUNCTION Tk(k,cosm)


    IMPLICIT NONE
    REAL :: Tk, k
    TYPE(cosmology) :: cosm

    !Transfer functions

    IF(cosm%itk==1) THEN
       Tk=1.
    ELSE IF(cosm%itk==2) THEN
       Tk=Tk_defw(k,cosm)
    ELSE IF(cosm%itk==3) THEN
       Tk=Tk_eh(k,cosm)
    ELSE IF(cosm%itk==4) THEN
       Tk=find_tk(k,cosm)
    END IF

  END FUNCTION Tk

  FUNCTION find_tk(k,cosm)


    IMPLICIT NONE
    REAL :: find_tk
    REAL :: kmin, kmax
    REAL, INTENT(IN) :: k
    INTEGER :: n
    TYPE(cosmology) :: cosm

    !Look-up and interpolation for T(k)

    n=SIZE(cosm%ktab)
    kmin=cosm%ktab(1)
    kmax=cosm%ktab(n)

    IF(k<kmin) THEN
       !For k<<keq Tk=1.
       find_tk=1.
    ELSE IF(k>kmax) THEN
       !Do some interpolation here based on knowledge of things at high k
       find_tk=cosm%tktab(n)*(log(k)/log(kmax))*((k/kmax)**(-2.))
    ELSE
       !Otherwise use the standard find algorithm
       find_tk=exp(find(log(k),log(cosm%ktab),log(cosm%tktab),3,3))
    END IF

  END FUNCTION find_tk

  FUNCTION find_pk(k,cosm)


    IMPLICIT NONE
    REAL :: find_pk
    REAL :: kmax
    REAL, INTENT(IN) :: k
    INTEGER :: n
    TYPE(cosmology) :: cosm

    !Look-up and interpolation for P(k)

    n=SIZE(cosm%ktab)
    kmax=cosm%ktab(n)

    IF(k>kmax) THEN
       !Do some interpolation here based on knowledge of things at high k
       find_pk=cosm%pktab(n)*((log(k)/log(kmax))**2.)*((k/kmax)**(cosm%n-1.))
    ELSE
       !Otherwise use the standard find algorithm
       find_pk=exp(find(log(k),log(cosm%ktab),log(cosm%pktab),3,3))
    END IF

  END FUNCTION find_pk

  FUNCTION Tk_eh(yy,cosm)

    ! the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
    ! remember I use k/h, whereas they use pure k, om_m is cdm + baryons


    IMPLICIT NONE

    REAL :: Tk_eh
    REAL, INTENT(IN) :: yy
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL*8 :: rk, e, thet, b1, b2, zd, ze, rd, re, rke, s, rks
    REAL*8 :: q
    REAL*8 :: y, g, ab
    REAL*8 :: a1, a2, ac
    REAL*8 :: bc
    REAL*8 :: f, fac
    REAL*8 :: c1, c2, tc
    REAL*8 :: bb, bn, ss, tb
    REAL*8 :: om_m, om_b, h

    om_m=cosm%om_m
    om_b=cosm%om_b
    h=cosm%h

    rk=yy*h

    e=exp(1.)

    thet=2.728/2.7
    b1=0.313*(om_m*h*h)**(-0.419)*(1+0.607*(om_m*h*h)**0.674)
    b2=0.238*(om_m*h*h)**0.223
    zd=1291.*(1+b1*(om_b*h*h)**b2)*(om_m*h*h)**0.251/(1.+0.659*(om_m*h*h)**0.828)
    ze=2.50e4*om_m*h*h/thet**4.
    rd=31500.*om_b*h*h/thet**4./zd
    re=31500.*om_b*h*h/thet**4./ze
    rke=7.46e-2*om_m*h*h/thet**2.
    s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1+sqrt(re)))
    rks=1.6*( (om_b*h*h)**0.52 ) * ( (om_m*h*h)**0.73 ) * (1.+(10.4*om_m*h*h)**(-0.95))

    q=rk/13.41/rke

    y=(1.+ze)/(1.+zd)
    g=y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)))
    ab=g*2.07*rke*s/(1.+rd)**(0.75)

    a1=(46.9*om_m*h*h)**0.670*(1+(32.1*om_m*h*h)**(-0.532))
    a2=(12.0*om_m*h*h)**0.424*(1+(45.0*om_m*h*h)**(-0.582))
    ac=(a1**(-om_b/om_m)) * (a2**(-(om_b/om_m)**3.))

    b1=0.944/(1+(458.*om_m*h*h)**(-0.708))
    b2=(0.395*om_m*h*h)**(-0.0266)
    bc=1./(1.+b1*((1.-om_b/om_m)**b2-1.))

    f=1./(1.+(rk*s/5.4)**4.)

    c1=14.2 + 386./(1.+69.9*q**1.08)
    c2=14.2/ac + 386./(1.+69.9*q**1.08)
    tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q)

    bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt((17.2*om_m*h*h)**2.+1.)
    bn=8.41*(om_m*h*h)**0.435
    ss=s/(1.+(bn/rk/s)**3.)**(1./3.)
    tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+(rk*s/5.2)**2.)
    IF((rk/rks**1.4)>7.) THEN
       fac=0.
    ELSE
       fac=exp(-(rk/rks)**1.4)
    END IF
    tb=(tb+ab*fac/(1.+(bb/rk/s)**3.))*sin(rk*ss)/rk/ss

    tk_eh=(om_b/om_m)*tb+(1-om_b/om_m)*tc

  END FUNCTION TK_EH

  FUNCTION Tk_DEFW(rk,cosm)


    IMPLICIT NONE
    REAL :: tk_DEFW, rk, gamma
    REAL :: rkeff, q, tk
    REAL*8 :: q8, tk8
    TYPE(cosmology) :: cosm

    !DEFW tramsfer function

    gamma=cosm%gamma

    rkeff=0.172+0.011*log(gamma/0.36)*log(gamma/0.36)
    q=1.e-20 + rk/gamma
    q8=1.e-20 + rkeff/gamma
    tk=1./(1.+(6.4*q+(3.0*q)**1.5+(1.7*q)**2)**1.13)**(1./1.13)
    tk8=1./(1.+(6.4*q8+(3.0*q8)**1.5+(1.7*q8)**2)**1.13)**(1./1.13)

    tk_defw=tk/tk8

  END FUNCTION Tk_DEFW

  FUNCTION p_lin(k,cosm)


    IMPLICIT NONE
    REAL :: p_lin, k
    TYPE(cosmology) :: cosm

    !This gives the linear power spectrum for the model in question
    !P(k) should have been previously normalised so as to get the amplitude 'A' correct

    IF(k==0.) THEN
       !If p_lin happens to be foolishly called for 0 mode (which should never happen, but might in integrals)
       p_lin=0.
    ELSE IF(k>1.e8) THEN
       !Avoids some issues if p_lin is called for very (absurdly) high k values
       !For some reason crashes can occur if this is the case
       p_lin=0.
    ELSE IF(cosm%itk==5) THEN
       !itk==5 means P(k) has been taken as an input file
       !In this case use the input P(k) file
       p_lin=(cosm%A**2.)*find_pk(k,cosm)
    ELSE
       !In this case look for the transfer function
       p_lin=(cosm%A**2.)*(Tk(k,cosm)**2.)*(k**(cosm%n+3.))
    END IF

  END FUNCTION p_lin

  FUNCTION p_2h(k,z,lut,plin,cosm)


    REAL :: p_2h
    REAL, INTENT(IN) :: k, plin, z
    REAL :: sigv, frac, q3, P13, P22, a, b, c, knl
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm

    !Produces the 'two-halo' power

    sigv=cosm%sigv
    frac=fdamp(z,cosm)
    IF(frac .NE. 0.) THEN
       p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2.)
    ELSE
       p_2h=plin
    END IF

  END FUNCTION p_2h

  FUNCTION p_1h(k,z,lut,cosm)


    IMPLICIT NONE
    REAL :: p_1h
    REAL, INTENT(IN) :: k, z
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: Dv, g, knl, fac, et, ks
    REAL, ALLOCATABLE :: integrand(:), wk(:)
    REAL :: sum
    INTEGER :: i
    REAL, PARAMETER :: pi=3.141592654

    !Does the one-halo power integral

    ALLOCATE(integrand(lut%n),wk(lut%n))
    integrand=0.
    wk=0.

    sum=0.

    !Only call eta once
    et=eta(z,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
       g=gnu(lut%nu(i))
       wk(i)=win(k*(lut%nu(i)**et),lut%rv(i),lut%c(i))
       integrand(i)=((lut%rv(i))**3.)*g*(wk(i)**2.)
    END DO

    !Carries out the integration
    sum=inttab(lut%nu,integrand,3)

    DEALLOCATE(integrand,wk)

    Dv=Delta_v(z,cosm)

    !These are just numerical factors from the 1-halo integral in terms of nu!
    p_1h=sum*2.*Dv*(k**3.)/(3.*pi)

    !Damping of the 1-halo term at very large scales
    ks=kstar(z,cosm)

    IF((k/ks)**2.>7.) THEN
       fac=0.
    ELSE
       fac=exp(-((k/ks)**2.))
    END IF

    p_1h=p_1h*(1.-fac)

  end function p_1h

  SUBROUTINE fill_sigtab(cosm)


    IMPLICIT NONE
    REAL :: rmin, rmax
    REAL, ALLOCATABLE :: rtab(:), sigtab(:)
    REAL :: r, sig
    INTEGER :: i, nsig
    TYPE(cosmology) :: cosm

    !This fills up tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigma(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is approximately power-law below and above these values of R   
    !This wouldn't be appropriate for models with a linear spectrum cut-off (e.g. WDM)

    !These must be not allocated before sigma calculations otherwise when sigma(r) is called
    !otherwise sigma(R) looks for the result in the tables
    IF(ALLOCATED(cosm%rtab)) DEALLOCATE(cosm%rtab)
    IF(ALLOCATED(cosm%sigtab)) DEALLOCATE(cosm%sigtab)   

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    rmin=1e-4
    rmax=1e3
    nsig=64

    IF(feedback) WRITE(*,*) 'SIGTAB: Filling sigma interpolation table'
    IF(feedback) WRITE(*,*) 'SIGTAB: Rmin:', rmin
    IF(feedback) WRITE(*,*) 'SIGTAB: Rmax:', rmax
    IF(feedback) WRITE(*,*) 'SIGTAB: Values:', nsig

    rmin=log(rmin)
    rmax=log(rmax)

    ALLOCATE(rtab(nsig),sigtab(nsig))

    DO i=1,nsig

       !Equally spaced r in log
       r=rmin+(rmax-rmin)*float(i-1)/float(nsig-1)
       r=exp(r)

       sig=sigma(r,cosm)

       rtab(i)=r
       sigtab(i)=sig

    END DO

    !Must be allocated after the sigtab calulation above
    ALLOCATE(cosm%rtab(nsig),cosm%sigtab(nsig))

    cosm%rtab=rtab
    cosm%sigtab=sigtab

    DEALLOCATE(rtab,sigtab)

    IF(feedback) WRITE(*,*) 'SIGTAB: Done'
    IF(feedback) WRITE(*,*)

  END SUBROUTINE fill_sigtab

  FUNCTION sigma(r,cosm)


    IMPLICIT NONE
    REAL :: sigma
    REAL, INTENT(IN) :: r
    REAL, ALLOCATABLE :: integrand(:)
    REAL :: rchange
    TYPE(cosmology), INTENT(IN) :: cosm
    
    !Sigma(R) calculation
    !Uses different intergation methods depending on r value

    rchange=1.e-2

    IF(ALLOCATED(cosm%sigtab) .EQV. .TRUE.) THEN

!       sigma=cosm%g*cosm%A*exp(find(log(r),log(cosm%rtab),log(cosm%sigtab),3,3))
       sigma=exp(find(log(r),log(cosm%rtab),log(cosm%sigtab),3,3))

    ELSE IF(r>=rchange) THEN

       !For large radii use the usual sigint
       sigma=sigint(r,cosm)

    ELSE IF(r<rchange) THEN
       
       !For smaller radii split the integral into two parts and sum the results
       sigma=sqrt(sigint1(r,cosm)+sigint2(r,cosm))

    END IF

  END FUNCTION sigma

  FUNCTION wk_tophat(x)

    IMPLICIT NONE
    REAL :: wk_tophat, x

    !The normlaised Fourier Transform of a top-hat
    !Taylor expansion used for low |x| to avoid cancellation problems

    IF(x<0.01) THEN
       wk_tophat=1.-(x**2.)/10.
    ELSE
       wk_tophat=3.*(sin(x)-x*cos(x))/(x**3.)
    END IF

  END FUNCTION wk_tophat

  FUNCTION inttab(x,y,iorder)

    IMPLICIT NONE
    REAL :: inttab
    REAL, INTENT(IN) :: x(:), y(:)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    REAL*8 :: sum
    INTEGER :: i, n, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    !Routine to integrate tables of data using the trapezium rule
    !Can either use linear, quadratic or cubic methods

    n=SIZE(x)

    IF(n .NE. SIZE(y)) STOP 'Tables must be of the same length'

    sum=0.d0

    IF(iorder==1) THEN

       !Sums over all Trapezia (a+b)*h/2
       DO i=1,n-1
          a=y(i+1)
          b=y(i)
          h=x(i+1)-x(i)
          sum=sum+(a+b)*h/2.d0
       END DO

    ELSE IF(iorder==2) THEN

       DO i=1,n-2

          x1=x(i)
          x2=x(i+1)
          x3=x(i+2)

          y1=y(i)
          y2=y(i+1)
          y3=y(i+2)

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
          q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
          q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

          !Takes value for first and last sections but averages over sections where you
          !have two independent estimates of the area
          IF(n==3) THEN
             sum=sum+q3-q1
          ELSE IF(i==1) THEN
             sum=sum+(q2-q1)+(q3-q2)/2.d0
          ELSE IF(i==n-2) THEN
             sum=sum+(q2-q1)/2.d0+(q3-q2)
          ELSE
             sum=sum+(q3-q1)/2.
          END IF

       END DO

    ELSE IF(iorder==3) THEN

       DO i=1,n-1

          !First choose the integers used for defining cubics for each section
          !First and last are different because the section does not lie in the *middle* of a cubic

          IF(i==1) THEN

             i1=1
             i2=2
             i3=3
             i4=4

          ELSE IF(i==n-1) THEN

             i1=n-3
             i2=n-2
             i3=n-1
             i4=n

          ELSE

             i1=i-1
             i2=i
             i3=i+1
             i4=i+2

          END IF

          x1=x(i1)
          x2=x(i2)
          x3=x(i3)
          x4=x(i4)

          y1=y(i1)
          y2=y(i2)
          y3=y(i3)
          y4=y(i4)

          CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

          !These are the limits of the particular section of integral
          xi=x(i)
          xf=x(i+1)

          qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
          qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

          sum=sum+qf-qi

       END DO

    END IF

    inttab=sum

  END FUNCTION inttab

  FUNCTION sigma_integrand(t,R,f,cosm)


    REAL :: sigma_integrand
    REAL, INTENT(IN) :: t, R
    REAL :: k, y, w_hat
    TYPE(cosmology), INTENT(IN) :: cosm

    INTERFACE
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    !Integrand to the sigma integral in terms of t. Defined by k=(1/t-1)/f(R) where f(R) is *any* function
    
    IF(t==0.) THEN
       !t=0 corresponds to k=infintiy when W(kR)=0.
       sigma_integrand=0.
    ELSE IF(t==1.) THEN
       !t=1 corresponds to k=0. when P(k)=0.
       sigma_integrand=0.
    ELSE
       !f(R) can be *any* function of R here to improve integration speed
       k=(-1.+1./t)/f(R)
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand=p_lin(k,cosm)*(w_hat**2.)/(t*(1.-t))
    END IF

  END FUNCTION sigma_integrand

  FUNCTION f_rapid(r)

    IMPLICIT NONE
    REAL :: f_rapid
    REAL, INTENT(IN) :: r
    REAL :: alpha

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error

    IF(r>1.e-2) THEN
       !alpha 0.3-0.5 works well
       alpha=0.5
    ELSE
       !If alpha=1 this goes tits up
       !alpha 0.7-0.9 works well
       alpha=0.8
    END IF

    f_rapid=r**alpha

  END FUNCTION f_rapid

  FUNCTION sigint(r,cosm)

    !Integrates between a and b until desired accuracy is reached!
    !Uses trapezium rule

    IMPLICIT NONE
    REAL :: r
    INTEGER :: i, j, jmax
    REAL :: sigint, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac
    REAL*8 :: sum1, sum2
    TYPE(cosmology), INTENT(IN) :: cosm

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    DO j=1,jmax

       n=ninit*2**(j-1)

       !Avoids the end-points where the integrand is 0 anyway
       DO i=2,n-1

          !x is defined on the interval 0 -> 1
          x=float(i-1)/float(n-1)
          sum2=sum2+sigma_integrand(x,r,f_rapid,cosm)

       END DO

       dx=1./float(n-1)
       sum2=sum2*dx
       sum2=sqrt(sum2)

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          sigint=sum2
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*)
          WRITE(*,*) 'SIGINT: r:', r
          WRITE(*,*) 'SIGINT: Integration timed out'
          WRITE(*,*)
          STOP
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION sigint

  FUNCTION sigint1(r,cosm)

    IMPLICIT NONE
    REAL :: r
    INTEGER :: i, j, jmax
    REAL :: sigint1, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac, xmin, xmax, k
    REAL*8 :: sum1, sum2
    TYPE(cosmology), INTENT(IN) :: cosm

    !Different sigma(R) integration routine

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    xmin=r/(r+r**.5)
    xmax=1.

    DO j=1,jmax

       n=ninit*2**(j-1)

       !Avoids the end-point where the integrand is 0 anyway
       DO i=1,n-1

          x=xmin+(xmax-xmin)*float(i-1)/float(n-1)

          IF(i==1 .OR. i==n) THEN
             fac=0.5
          ELSE
             fac=1.
          END IF

          k=(-1.+1./x)/r**.5
          sum2=sum2+fac*p_lin(k,cosm)*(wk_tophat(k*r)**2.)/(x*(1.-x))

       END DO

       dx=(xmax-xmin)/float(n-1)
       sum2=sum2*dx

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          sigint1=sum2
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*)
          WRITE(*,*) 'SIGINT1: r:', r
          WRITE(*,*) 'SIGINT1: Integration timed out'
          WRITE(*,*)
          STOP
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION sigint1

  FUNCTION sigint2(r,cosm)

    IMPLICIT NONE
    REAL :: r
    INTEGER :: i, j, jmax
    REAL :: sigint2, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac, xmin, xmax, A
    REAL*8 :: sum1, sum2
    TYPE(cosmology), INTENT(IN) :: cosm
    
    !Another different sigma(R) integration routine

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    !How far to go out in 1/r units for integral
    A=10.

    xmin=1./r
    xmax=A/r

    DO j=1,jmax

       n=ninit*2**(j-1)

       DO i=1,n

          x=xmin+(xmax-xmin)*float(i-1)/float(n-1)

          IF(i==1 .OR. i==n) THEN
             fac=0.5
          ELSE
             fac=1.
          END IF

          !Integrate linearly in k for the rapidly oscillating part
          sum2=sum2+fac*p_lin(x,cosm)*(wk_tophat(x*r)**2.)/x

       END DO

       dx=(xmax-xmin)/float(n-1)
       sum2=sum2*dx

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          sigint2=sum2
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*)
          WRITE(*,*) 'SIGINT2: r:', r
          WRITE(*,*) 'SIGINT2: Integration timed out'
          WRITE(*,*)
          STOP
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION sigint2

  FUNCTION win(k,rv,c)

    IMPLICIT NONE
    REAL :: win, k, rv, c

    !Calls the analytic Fourier Transform of the NFW profile

    win=winnfw(k,rv,c)

    !Correct for the case of disasters (a bit sloppy, not sure if this is ever used)
    IF(win>1.) win=1.
    IF(win<0.) win=0.

  END FUNCTION win

  FUNCTION winnfw(k,rv,c)

    IMPLICIT NONE
    REAL :: winnfw
    REAL, INTENT(IN) :: k, rv, c
    REAL :: si1, si2, ci1, ci2, ks
    REAL :: p1, p2, p3

    !The analytic Fourier Transform of the NFW

    ks=k*rv/c

    si1=si(ks)
    si2=si((1.+c)*ks)
    ci1=ci(ks)
    ci2=ci((1.+c)*ks)

    p1=cos(ks)*(ci2-ci1)
    p2=sin(ks)*(si2-si1)
    p3=sin(ks*c)/(ks*(1.+c))

    winnfw=p1+p2-p3
    winnfw=winnfw/mass(c)

  END FUNCTION winnfw

  FUNCTION mass(c)

    !This calculates the (normalised) mass of a halo of concentration c
    !This is mass with some factors missing (4*pi, rs ?)

    IMPLICIT NONE
    REAL :: mass, c

    mass=log(1.+c)-c/(1.+c)

  END FUNCTION mass

  FUNCTION gnu(nu)

    IMPLICIT NONE
    REAL :: gnu, nu

    !Mass function

    gnu=gst(nu)

  END FUNCTION gnu

  FUNCTION gst(nu)

    !Sheth Tormen mass function!

    IMPLICIT NONE
    REAL :: nu, gst
    REAL :: p, q

    p=0.3
    q=0.707

    gst=0.21616*(1.+((q*nu*nu)**(-p)))*exp(-q*nu*nu/2.)

  END FUNCTION gst

  FUNCTION hubble2(z,cosm)

    !This calculates the dimensionless squared hubble parameter at redshift z!
    !and it ignores contributions from radiation (not accurate at high z, but consistent with simulations)!


    IMPLICIT NONE
    REAL :: hubble2, z
    REAL :: om_m, om_v, w
    TYPE(cosmology) :: cosm

    om_m=cosm%om_m
    om_v=cosm%om_v
    w=cosm%w

    hubble2=(om_m*(1.+z)**3.)+om_v*((1.+z)**(3.*(1.+w)))+((1.-om_m-om_v)*(1.+z)**2.)

  END FUNCTION hubble2

  FUNCTION omega_m(z,cosm)

    !This calculates Omega_m variations with z!


    IMPLICIT NONE
    REAL :: omega_m, z
    REAL :: om_m
    TYPE(cosmology) :: cosm

    om_m=cosm%om_m

    omega_m=(om_m*(1.+z)**3.)/hubble2(z,cosm)

  END FUNCTION omega_m

  FUNCTION omega_v(z,cosm)

    !This calculates Omega_v variations with z for any w


    IMPLICIT NONE
    REAL :: omega_v, z
    REAL :: om_v, w
    TYPE(cosmology) :: cosm

    om_v=cosm%om_v
    w=cosm%w

    omega_v=om_v*((1.+z)**(3.*(1.+w)))/hubble2(z,cosm)

  END FUNCTION omega_v

  FUNCTION grow(z,cosm)


    IMPLICIT NONE
    REAL :: grow, z, a
    TYPE(cosmology) :: cosm

    !Computes the growth function

    IF(z==0.) THEN
       grow=1.
    ELSE
       a=1./(1.+z)
       grow=grow_int(1.,a,0.001,cosm)
    END IF

  END FUNCTION grow

  FUNCTION grow_int(a,b,acc,cosm)


    IMPLICIT NONE
    INTEGER :: i, j, jmax
    REAL :: grow_int, a, b, acc, dx
    INTEGER :: nint
    REAL :: x, fac, func, this, gam
    REAL*8 :: sum1, sum2
    TYPE(cosmology) :: cosm

    !Integrates the Linder (2005) expression to find the growth function

    sum1=0.d0
    sum2=0.d0

    jmax=20

    DO j=1,jmax

       nint=10.*(2.**j)

       DO i=1,nint

          x=a+(b-a)*((float(i)-1)/(float(nint)-1))

          IF(i==1 .OR. i==nint) THEN
             !multiple of 1 for beginning and end and multiple of 2 for middle points!
             fac=1.
          ELSE
             fac=2.
          END IF

          IF(cosm%w<-1.) THEN
             gam=0.55+0.02*(1.+cosm%w)
          ELSE IF(cosm%w>-1) THEN
             gam=0.55+0.05*(1.+cosm%w)
          ELSE
             gam=0.55
          END IF

          func=(omega_m(-1.+1./x,cosm)**gam)/x

          sum2=sum2+fac*func

       END DO

       dx=((b-a)/(float(nint)-1.))
       sum2=sum2*dx/2.

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          grow_int=exp(sum2)
          EXIT
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION grow_int

  FUNCTION dispint(cosm)


    IMPLICIT NONE
    REAL :: dispint
    REAL*8 :: sum
    REAL :: dtheta, k, theta, oldsum, acc
    REAL, PARAMETER :: pi=3.141592654
    INTEGER :: i, j, n, ninit, jmax
    TYPE(cosmology), INTENT(IN) :: cosm

    !Integration routine to calculate sigma_v

    acc=0.001
    ninit=100
    jmax=30

    DO j=1,jmax

       n=ninit*(2**(j-1))

       sum=0.d0
       dtheta=1./float(n)

       DO i=2,n-1

          !theta converts integrand to 0->1 range
          !Values at the end points are 0 so removed for convenience
          theta=float(i-1)/float(n-1)
          k=(-1.+1./theta)          
          sum=sum+((1.+k)**2.)*p_lin(k,cosm)/(k**3.)

       END DO

       sum=sum*dtheta

       IF(j>1 .AND. ABS(-1.+sum/oldsum)<acc) THEN  
          dispint=sum/3.
          EXIT
       ELSE
          oldsum=sum
       END IF

    END DO

  END FUNCTION dispint

  FUNCTION si(x)

    IMPLICIT NONE
    REAL :: si, x
    REAL*8 :: x2, y, f, g, si8
    REAL*8, PARAMETER :: pi=3.1415926535897932384626433d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.

    IF(ABS(x)<=4.) THEN

       x2=x*x

       si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3&
            +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10&
            +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
            (1.+x2*(1.01162145739225565d-2 +x2*(4.99175116169755106d-5+&
            x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13&
            +x2*(3.21107051193712168d-16)))))))

       si=si8

    ELSE IF(ABS(x)>4.) THEN

       y=1.d0/(x*x)

       f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 +&
            y*(2.37750310125431834034d7 +y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10 &
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + &
            y*(1.00795182980368574617d13 + y*(4.94816688199951963482d12 +&
            y*(-4.94701168645415959931d11)))))))))))/ (x*(1. +y*(7.46437068161927678031d2 +&
            y*(1.97865247031583951450d5 +y*(2.41535670165126845144d7 + &
            y*(1.47478952192985464958d9 + y*(4.58595115847765779830d10 +&
            y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 + &
            y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))


       g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + &
            y*(3.12557570795778731d7 + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 +&
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 +y*(1.81004487464664575d13 +&
            y*(6.43291613143049485d12 +y*(-1.36517137670871689d12)))))))))))/&
            (1. + y*(8.19595201151451564d2 +y*(2.40036752835578777d5 + y*(3.26026661647090822d7 &
            + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 + y*(1.39866710696414565d12 &
            + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 +y*(3.99653257887490811d13))))))))))

       si=pi/2.d0-f*cos(x)-g*sin(x)

    END IF

  END FUNCTION si

  FUNCTION ci(x)

    IMPLICIT NONE
    REAL :: ci, x
    REAL*8 :: x2, y, f, g, ci8
    REAL*8, PARAMETER :: em_const=0.577215664901532861d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.

    IF(ABS(x)<=4.) THEN

       x2=x*x

       ci8=em_const+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4&
            +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11&
            +x2*(-9.93728488857585407d-15)))))))/ (1.+x2*(1.1592605689110735d-2+&
            x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+&
            x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

       ci=ci8

    ELSE IF(ABS(x)>4.) THEN

       y=1./(x*x) 

       f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 + &
            y*(2.37750310125431834034d7 +y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10&
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + y*(1.00795182980368574617d13&
            + y*(4.94816688199951963482d12 +y*(-4.94701168645415959931d11)))))))))))/&
            (x*(1. +y*(7.46437068161927678031d2 +y*(1.97865247031583951450d5 +&
            y*(2.41535670165126845144d7 + y*(1.47478952192985464958d9 + &
            y*(4.58595115847765779830d10 +y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 &
            + y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))   

       g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + y*(3.12557570795778731d7&
            + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 +&
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 +&
            y*(1.81004487464664575d13 + y*(6.43291613143049485d12 +y*(-1.36517137670871689d12)))))))))))&
            / (1. + y*(8.19595201151451564d2 +y*(2.40036752835578777d5 +&
            y*(3.26026661647090822d7 + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 &
            + y*(1.39866710696414565d12 + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 +y*(3.99653257887490811d13))))))))))

       ci=f*sin(x)-g*cos(x)

    END IF

  END FUNCTION ci

  FUNCTION derivative_table(x,xin,yin,iorder,imeth)

    IMPLICIT NONE
    REAL :: derivative_table
    REAL, INTENT(IN) :: x, xin(:), yin(:)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !Finds the derivative f'(x) given tables x, f(x)

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !imeth = 1 => find x in xtab by crudely searching
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    n=SIZE(xtab)

    maxorder=3
    maxmethod=3

    n=SIZE(xin)
    IF(n .NE. SIZE(yin)) STOP 'FIND: Tables not of the same size'
    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab)
       CALL reverse(ytab)
    END IF

    IF(iorder<1) STOP 'FIND: find order not specified correctly'
    IF(iorder>maxorder) STOP 'FIND: find order not specified correctly'
    IF(imeth<1) STOP 'FIND: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) STOP 'FIND: Method of finding within a table not specified correctly'

    IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x2=xtab(2)
          x1=xtab(1)

          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-1)) THEN

          x2=xtab(n)
          x1=xtab(n-1)

          y2=ytab(n)
          y1=ytab(n-1)

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x2=xtab(i+1)
          x1=xtab(i)

          y2=ytab(i+1)
          y1=ytab(i)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       derivative_table=a

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND_QUADRATIC: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x3=xtab(3)
             x2=xtab(2)
             x1=xtab(1)

             y3=ytab(3)
             y2=ytab(2)
             y1=ytab(1)

          ELSE IF (x>=xtab(n-1)) THEN

             x3=xtab(n)
             x2=xtab(n-1)
             x1=xtab(n-2)

             y3=ytab(n)
             y2=ytab(n-1)
             y1=ytab(n-2)

          END IF

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          derivative_table=2.*a*x+b

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          !In this case take the average of two separate quadratic spline values

          derivative_table=0.

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          derivative_table=derivative_table+(2.*a*x+b)/2.

          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          derivative_table=derivative_table+(2.*a*x+b)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND_CUBIC: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x4=xtab(4)
          x3=xtab(3)
          x2=xtab(2)
          x1=xtab(1)

          y4=ytab(4)
          y3=ytab(3)
          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-2)) THEN

          x4=xtab(n)
          x3=xtab(n-1)
          x2=xtab(n-2)
          x1=xtab(n-3)

          y4=ytab(n)
          y3=ytab(n-1)
          y2=ytab(n-2)
          y1=ytab(n-3)

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       derivative_table=3.*a*(x**2.)+2.*b*x+c

    END IF

  END FUNCTION derivative_table

  FUNCTION find(x,xin,yin,iorder,imeth)

    IMPLICIT NONE
    REAL :: find
    REAL, INTENT(IN) :: x, xin(:), yin(:)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !Interpolation routine.

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    maxorder=3
    maxmethod=3

    n=SIZE(xin)
    IF(n .NE. SIZE(yin)) STOP 'FIND: Tables not of the same size'
    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab)
       CALL reverse(ytab)
    END IF

    IF(iorder<1) STOP 'FIND: find order not specified correctly'
    IF(iorder>maxorder) STOP 'FIND: find order not specified correctly'
    IF(imeth<1) STOP 'FIND: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) STOP 'FIND: Method of finding within a table not specified correctly'

    IF(x<xtab(1)) THEN

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(x>xtab(n)) THEN

       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x1=xtab(1)
          x2=xtab(2)

          y1=ytab(1)
          y2=ytab(2)

       ELSE IF (x>=xtab(n-1)) THEN

          x1=xtab(n-1)
          x2=xtab(n)

          y1=ytab(n-1)
          y2=ytab(n)

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i)
          x2=xtab(i+1)

          y1=ytab(i)
          y2=ytab(i+1)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x1=xtab(1)
             x2=xtab(2)
             x3=xtab(3)

             y1=ytab(1)
             y2=ytab(2)
             y3=ytab(3)

          ELSE IF (x>=xtab(n-1)) THEN

             x1=xtab(n-2)
             x2=xtab(n-1)
             x3=xtab(n)

             y1=ytab(n-2)
             y2=ytab(n-1)
             y3=ytab(n)

          END IF

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          find=a*(x**2.)+b*x+c

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          !In this case take the average of two separate quadratic spline values

          find=0.

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          find=find+(a*(x**2.)+b*x+c)/2.

          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          find=find+(a*(x**2.)+b*x+c)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x1=xtab(1)
          x2=xtab(2)
          x3=xtab(3)
          x4=xtab(4)        

          y1=ytab(1)
          y2=ytab(2)
          y3=ytab(3)
          y4=ytab(4)

       ELSE IF (x>=xtab(n-2)) THEN

          x1=xtab(n-3)
          x2=xtab(n-2)
          x3=xtab(n-1)
          x4=xtab(n)

          y1=ytab(n-3)
          y2=ytab(n-2)
          y3=ytab(n-1)
          y4=ytab(n)

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       find=a*x**3.+b*x**2.+c*x+d

    END IF

  END FUNCTION find

  FUNCTION linear_table_integer(x,xtab)

    IMPLICIT NONE
    INTEGER :: linear_table_integer
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: n
    REAL :: x1, x2, xn
    REAL :: acc

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
    !Assumes table is organised linearly (care for logs)

    n=SIZE(xtab)
    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)

    !Test for linear table
    acc=0.001

    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

  END FUNCTION linear_table_integer

  FUNCTION search_int(x,xtab)

    IMPLICIT NONE
    INTEGER :: search_int
    INTEGER :: i, n
    REAL, INTENT(IN) :: x, xtab(:)

    !Searches for the point in the table brute force.
    !This is usually a stupid thing to do

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
       IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

  END FUNCTION search_int

  FUNCTION int_split(x,xtab)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: i1, i2, imid, n
    INTEGER :: int_split

    !Finds the position of the value in the table by continually splitting it in half

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) STOP 'INT_SPLIT: table in wrong order'

    i1=1
    i2=n

    DO

       imid=NINT((i1+i2)/2.)

       IF(x<xtab(imid)) THEN
          i2=imid
       ELSE
          i1=imid
       END IF

       IF(i2==i1+1) EXIT

    END DO

    int_split=i1

  END FUNCTION int_split

  SUBROUTINE fit_line(a1,a0,x1,y1,x2,y2)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2

    !Given xi, yi i=1,2 fits a line between these points

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

  END SUBROUTINE fit_line

  SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3

    !Given xi, yi i=1,2,3 fits a quadratic between these points

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.)-a1*x1

  END SUBROUTINE fit_quadratic

  SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3

    !Given xi, yi i=1,2,3,4 fits a cubic between these points

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2.+x4*x1+x1**2.)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3.-b*x1**2.-c*x1

  END SUBROUTINE fit_cubic

  SUBROUTINE reverse(arry)

    IMPLICIT NONE
    INTEGER :: n, i
    REAL, ALLOCATABLE :: hold(:)
    REAL :: arry(:)

    !This reverses the contents of 'arry'

    n=SIZE(arry)

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

  END SUBROUTINE reverse

  FUNCTION file_length(file_name)

    IMPLICIT NONE
    CHARACTER(len=64) :: file_name
    INTEGER ::n, file_length
    REAL :: data

    !Finds the length of a file

    OPEN(7,file=file_name)
    n=0
    DO
       n=n+1
       READ(7,*, end=301) data
    END DO

    !301 is just the label to jump to when the end of the file is reached

301 CLOSE(7)  

    file_length=n-1

  END FUNCTION file_length

END MODULE MHM


