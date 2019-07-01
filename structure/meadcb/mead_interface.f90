module mead_settings_mod
    type mead_settings
        logical :: noisy
        logical :: feedback
        logical :: one_baryon_param
        !real(8) :: kmin, kmax
        !integer :: nk

        real(8) :: zmin, zmax
        integer :: nz

        character(len=256) :: linear_ps_section_name, output_section_name

    end type mead_settings

end module mead_settings_mod

function setup(options) result(result)
    use mead_settings_mod
    use cosmosis_modules
    implicit none

    integer(cosmosis_block), value :: options
    integer(cosmosis_status) :: status
    type(mead_settings), pointer :: settings
    type(c_ptr) :: result
    status = 0
    
    allocate(settings)
    status = status + datablock_get_logical_default(options, option_section, "feedback", .false., settings%feedback)

    status = status + datablock_get_double_default(options, option_section, "zmin", real(0.0, kind=8), settings%zmin)
    status = status + datablock_get_double_default(options, option_section, "zmax", real(3.0, kind=8), settings%zmax)
    status = status + datablock_get_int_default(options, option_section, "nz", -1, settings%nz)

    status = status + datablock_get_logical_default(options, option_section, "one_baryon_parameter", .false., settings%one_baryon_param)

    status = datablock_get_string_default(options, option_section, "input_section_name", matter_power_lin_section, settings%linear_ps_section_name)
    status = datablock_get_string_default(options, option_section, "output_section_name", matter_power_nl_section, settings%output_section_name)


    if (status .ne. 0) then
        write(*,*) "One or more parameters not found for hmcode"
        stop
    endif

    result = c_loc(settings)

end function setup


function execute(block,config) result(status)
    use mead_settings_mod
    use cosmosis_modules
    use mhmcamb

    implicit none

    logical :: feedback
    integer(cosmosis_block), value :: block
    integer(cosmosis_status) :: status
    type(c_ptr), value :: config
    type(mead_settings), pointer :: settings    
    integer, parameter :: LINEAR_SPACING = 0
    integer, parameter :: LOG_SPACING = 1
    character(*), parameter :: cosmo = cosmological_parameters_section
    character(*), parameter :: halo = halo_model_parameters_section

    real(4) :: p1h, p2h,pfull, plin, z
    integer :: i,j,nk,nz, massive_nu
    REAL, ALLOCATABLE :: k(:), ztab(:)
    TYPE(HM_cosmology) :: cosi
    TYPE(HM_tables) :: lut
    !CosmoSIS supplies double precision - need to convert
    real(8) :: om_m, om_v, om_b, h, w, n_s, om_nu
    real(8) :: wa, t_cmb
    real(8), ALLOCATABLE :: k_in(:), z_in(:), p_in(:,:)
    real(8), ALLOCATABLE :: k_out(:), z_out(:), p_out(:,:)
    real(8) :: halo_as, halo_eta0

    HM_verbose = .False.
    imead = 1
    status = 0
    call c_f_pointer(config, settings)

    feedback = settings%feedback

    !Fill in the cosmology parameters. We need to convert from CosmoSIS 8-byte reals
    !to HMcode 4-byte reals, hence the extra bit
    status = status + datablock_get(block, cosmo, "omega_m", om_m)
    status = status + datablock_get(block, cosmo, "omega_lambda", om_v)
    status = status + datablock_get(block, cosmo, "omega_b", om_b)
    status = status + datablock_get_double_default(block, cosmo, "omega_nu", 0.0D0, om_nu)
    status = status + datablock_get_double_default(block, cosmo, "t_cmb", 2.726D0, t_cmb)
    status = status + datablock_get_int_default(block, cosmo, "massive_nu", 3, massive_nu)
    status = status + datablock_get(block, cosmo, "h0", h)
    status = status + datablock_get(block, cosmo, "n_s", n_s)
    status = status + datablock_get_double_default(block, cosmo, "w", -1.0D0, w)
    status = status + datablock_get_double_default(block, cosmo, "wa", 0.0D0, wa)
    status = status + datablock_get_double_default(block, halo, "A", 3.13D0, halo_as)
    
    if(.not. settings%one_baryon_param) then
        status = status + datablock_get_double_default(block, halo, "eta_0", 0.603D0, halo_eta0)
    else
        halo_eta0 = 1.03-0.11*halo_as
    endif

    if (status .ne. 0 ) then
        write(*,*) "Error reading parameters for Mead code"
        return
    endif

 
    status = status + datablock_get_double_grid(block, settings%linear_ps_section_name, &
        "k_h", k_in, "z", z_in, "p_k", p_in)

    if (status .ne. 0 ) then
        write(*,*) "Error reading P(k,z) for Mead code"
        return
    endif

    ! doing the job of initialise_HM_cosmology
    ! 1. get linear pk(z) (for each z)
    ! 2. fill sigma table (for each z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --cosi ->cosm, As->Abary, ktab->k_plin(?), eta_0->eta0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --hmcode only reads pk_lin at z = 0.
    !when doing mead, camb needs to start from 0.0<- oldversion
    !new hmcode take pkz input to label redshift of p_lin. 

    !Find the index of z where z==0

    !Copy in P(k) from the right part of P(k,z)
    nk = size(k_in)

    if(settings%nz > 0) then
        nz = settings%nz
    else
        nz = size(z_in)
    endif
 
    ! doing the job of assign_HM_cosmology
    cosi%om_m=om_m !The halo modelling in this version knows about neutrino
    cosi%om_v=om_v
    cosi%f_nu=om_nu/cosi%om_m
    cosi%Tcmb=t_cmb
    cosi%h=h
    cosi%w=w
    cosi%ns=n_s
    cosi%wa = wa
    cosi%Nnu=massive_nu
    !testing massive_nu effect

    cosi%eta_baryon = halo_eta0
    cosi%A_baryon = halo_as
    cosi%nk =nk

    !Fill growth function table (only needs to be done once)
    CALL fill_growtab(cosi)
    !getting linear growth
    
    !And get the cosmo power spectrum, again as double precision
    !Also the P is 2D as we get z also

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(k(nk))
    ALLOCATE(ztab(nz))
        
    !Set the output ranges in k and z
    !this was not here previously, take care of output k and pk
    k = k_in
    if(settings%nz > 0) then
        CALL fill_table(real(settings%zmin),real(settings%zmax),ztab,settings%nz)
    else
        ztab = z_in
    endif
        
    !Fill table for output power
    ALLOCATE(p_out(nk,nz))
 
    !Loop over redshifts
    DO j=1,nz
        CALL fill_plintab_cosmosis(j,cosi,real(k_in),real(z_in),real(p_in),size(k_in),size(z_in))
        CALL fill_sigtab(cosi)
        !Sets the redshift
        z=ztab(j)

        !Initialisation for the halomodel calcualtion
        !Also normalises power spectrum (via sigma_8)
        !and fills sigma(R) tables 
        CALL halomod_init(z,lut,cosi)

        !Loop over k values
        DO i=1,nk
            plin = p_in(i,j)*(k(i)**3.0) / (2.*(pi**2.)) 
            !plin=p_lin(k(i),z,0,cosi) using function in hmcode to calculate p_lin.
            !not necessary any more, p_lin result and pk_lin in cosmosis has maximum diff ~1e-6
            CALL halomod(k(i),z,p1h,p2h,pfull,plin,lut,cosi)
            !This outputs k^3 P(k).  We convert back. (need to check for new version) 
            ! note, i,j are interchanged with respect to hm main program
            p_out(i,j)=pfull / (k(i)**3.0) * (2.*(pi**2.)) 
        END DO

        IF(j==1) THEN
            if (settings%feedback) WRITE(*,fmt='(A5,A7)') 'i', 'z'
            if (settings%feedback) WRITE(*,fmt='(A13)') '   ============'
        END IF
        if (settings%feedback) WRITE(*,fmt='(I5,F8.3)') j, ztab(j)
    END DO

    !convert to double precision
    allocate(k_out(nk))
    allocate(z_out(nz))
    k_out = k
    z_out = ztab
    !Convert k to k/h to match other modules
    !Output results to cosmosis
    status = datablock_put_double_grid(block, settings%output_section_name, "k_h", k_out, "z", z_out, "p_k", p_out)

    !Free memory
    deallocate(k)
    deallocate(ztab)
    deallocate(p_out)
    deallocate(k_in)
    deallocate(z_in)
    deallocate(p_in)
    deallocate(k_out)
    deallocate(z_out)
    call deallocate_LUT(lut)
    IF(ALLOCATED(cosi%k_plin)) DEALLOCATE(cosi%k_plin)
    IF(ALLOCATED(cosi%plin)) DEALLOCATE(cosi%plin)
    IF(ALLOCATED(cosi%plinc)) DEALLOCATE(cosi%plinc)   
    IF(ALLOCATED(cosi%r_sigma)) DEALLOCATE(cosi%r_sigma)
    IF(ALLOCATED(cosi%sigma)) DEALLOCATE(cosi%sigma)
    IF(ALLOCATED(cosi%a_growth)) DEALLOCATE(cosi%a_growth)
    IF(ALLOCATED(cosi%growth)) DEALLOCATE(cosi%growth)
    
end function execute
