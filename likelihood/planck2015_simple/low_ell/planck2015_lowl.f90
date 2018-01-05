module planck2015_lowl_tools
    use comm_gauss_br_mod, only : comm_gauss_br_initialize_object, &
    comm_gauss_br_deallocate_object, comm_gauss_br_compute_lnL, comm_gauss_br
    implicit none

    type planck2015_low_settings
        integer handle
        integer lmax
    end type

contains


end module planck2015_lowl_tools


function setup(options) result(result)
    use cosmosis_modules
    use planck2015_lowl_tools
    implicit none
    integer(cosmosis_block), value :: options
    integer(cosmosis_status) :: status
    type(planck2015_low_settings), pointer :: settings
    type(c_ptr) :: result
    integer :: lmin, lmax, delta_l
    integer, parameter :: default_lmin=2, default_lmax=29, default_delta_l=29, default_handle=1
    character(512) :: gaussfile
    character(512) :: default_gaussfile
    character(512) :: cosmosis_dir

    !Make space for the configuration
    allocate(settings)
    status = 0
    status = status + datablock_get_int_default(options, option_section, "handle", default_handle, settings%handle)
    status = status + datablock_get_int_default(options, option_section, "lmin", default_lmin, lmin)
    status = status + datablock_get_int_default(options, option_section, "lmax", default_lmax, settings%lmax)
    status = status + datablock_get_int_default(options, option_section, "delta_l", default_delta_l, delta_l)

    if (status .ne. 0) then
        write(*,*) "Error setting handle, lmin, lmax, or delta_l - these should all be integers or left unset"
        stop
    endif

    CALL get_environment_variable("COSMOSIS_SRC_DIR", cosmosis_dir)
    default_gaussfile = trim(cosmosis_dir) // "/cosmosis-standard-library/likelihood/planck2015/data/commander_rc2_v1.1_l2_29_B.clik/clik/lkl_0/_external/sigma.fits" 

    gaussfile = ""
    status = status + datablock_get_string_default(options, option_section, "gaussfile", default_gaussfile, gaussfile)
    if (status .ne. 0) then
        write(*,*) "Error setting gaussfile parameter - should be a string"
        stop
    endif

    write(*,*) "Loading Planck low-ell data from: ", trim(gaussfile)

    call comm_gauss_br_initialize_object(gaussfile, lmin, settings%lmax, delta_l, settings%handle)

    result = c_loc(settings)
end function setup


function execute(block, config) result(status)
    use cosmosis_modules
    use planck2015_lowl_tools

    implicit none
    integer(cosmosis_block), value :: block
    integer(cosmosis_status) :: status
    type(c_ptr), value :: config
    type(planck2015_low_settings), pointer :: settings
    real(8) :: like
    real(8), dimension(:), allocatable :: tt_in, tt
    real(8), dimension(:), allocatable :: x
    integer, dimension(:), allocatable :: ell
    real(8) :: A_planck
    integer :: n_ell, l


    call c_f_pointer(config, settings)
    status = 0
    status = datablock_get_double_array_1d(block, cmb_cl_section, "tt", tt_in, n_ell)
    status = datablock_get_int_array_1d(block, cmb_cl_section, "ell", ell, n_ell)

    if (status .ne. 0) then
        write(*,*) "tt or ell not calculated for quick_commander"
        status = 1
        return
    endif

    if (.not. (ell(n_ell) .gt. settings%lmax)) then
        deallocate(tt_in, ell)
        write(*,*) "ell not calculated high enough for the "
        status = 2
        return
    endif

    status = datablock_get_double(block, "planck", "A_planck", A_planck)
    if (status .ne. 0) then
        deallocate(tt_in, ell)
        write(*,*) "Calibration parameter A_planck not found in the [planck] section"
        status = 1
        return
    endif


    !Read cls
    allocate(tt(2:settings%lmax))
    do l=2,settings%lmax
        tt(l) = tt_in(l-ell(1)+1) / A_planck**2
    enddo


    like = comm_gauss_br_compute_lnL(tt, settings%handle, x)
  
    deallocate(tt_in, ell, tt)
 
    status = status + datablock_put_double(block, likelihoods_section, &
        "planck2015_lowl_like", like)

    status = status + datablock_put_double_array_1d(block, data_vector_section,&
                     "planck2015_lowl_theory", x)
    status = status + datablock_put_double_array_2d(block, data_vector_section,&
                     "planck2015_lowl_inverse_covariance", comm_gauss_br(settings%handle)%cov)

    write(*,*) status
end function execute



function cleanup(config) result(status)
    use cosmosis_modules
    use planck2015_lowl_tools
    type(c_ptr), value :: config
    type(planck2015_low_settings), pointer :: settings  
    integer(cosmosis_status) :: status

    !Free memory allocated in the setup function
    call c_f_pointer(config, settings)
    call comm_gauss_br_deallocate_object(settings%handle)
    deallocate(settings)



    status = 0

end function cleanup
