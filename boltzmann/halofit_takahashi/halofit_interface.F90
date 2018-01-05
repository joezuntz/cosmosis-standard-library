module halofit_interface_tools
use cosmosis_modules
use transfer
use nonlinear
implicit none

	type halofit_settings
		real(8) :: kmin, kmax
		integer :: nk
	end type


contains


function load_matter_power(block, PK) result(status)
	use cosmosis_modules
	use nonlinear
	integer(cosmosis_block) :: block
	integer(cosmosis_status) :: status
	type(MatterPowerData) :: PK
	real(dl), allocatable, dimension(:) :: k, z
	real(dl), allocatable, dimension(:,:) :: P

	logical :: k_changes_fastest
	integer nkz
	integer i,j
	
	!Get the data columns from the fits data
	status = 0

	!Load k, z, P
	status = datablock_get_double_grid(block, matter_power_lin_section, &
		"K_H", k, "Z", z, "P_K", P)

	if (status .ne. 0) then
		write(*,*) "Could not find K_H, Z, or P_K in block"
		return
	endif

	!Fill in Halofit's data structure
	PK%num_k = size(k)
	PK%num_z = size(z)
	call allocate_matterpower(PK)
	PK%log_kh = log(k)
	PK%redshifts = z
	PK%matpower = log(P)

	!Allocate memory
	call MatterPowerdata_getsplines(PK)

	!Clean up
	deallocate(k, z, P)
	
end function

end module halofit_interface_tools

function setup(options) result(result)
	use cosmosis_modules
	use halofit_interface_tools
	implicit none
	integer(cosmosis_block), value :: options
	integer(cosmosis_status) :: status
	type(halofit_settings), pointer :: settings
	type(c_ptr) :: result
	allocate(settings)

	!Set default kmin and kmax to -1. We'll later set to the k_min and
	!k_max of the linear P(k), if not specified here
	status = 0
	status = status + datablock_get_double_default(options, option_section, "kmin", -1.0d0, settings%kmin)
	status = status + datablock_get_double_default(options, option_section, "kmax", -1.0d0, settings%kmax)
	status = status + datablock_get_int_default(options, option_section, "nk", 200, settings%nk)

	if (status .ne. 0) then 
		write(*,*) "Failed setup of halofit!", status
		stop
	endif
	result = c_loc(settings)
end function setup

function execute(block, config) result(status)
	use halofit_interface_tools
	use cosmosis_modules
	use lambdageneral
	implicit none
	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	type(c_ptr), value :: config
	type(halofit_settings), pointer :: settings	
	type(MatterPowerdata) :: PK
	type(MatterPowerdata) :: PK_NL
	real(dl), dimension(:,:), allocatable :: p
	integer ik, nk
	real(dl) :: kmin, kmax, zmin, zmax, log_kmin, log_kmax
	real(dl) :: omega_baryon, omega_matter, omega_nu
	real(dl), allocatable, dimension(:) :: k
	integer iz

	status = 0
	call c_f_pointer(config, settings)

	
	!Set Halofit internal numbers
	status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_B", omega_baryon)
	status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_M", omega_matter)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "OMEGA_NU", 0.0D0, omega_nu)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "W", -1.0D0, w_lam)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "WA", 0.0D0, wa_ppf)
	cp%omegac = omega_matter - omega_baryon - omega_nu
	cp%omegab = omega_baryon
	cp%omegan = omega_nu
	cp%omegav = 1 - omega_matter

    if (status .ne. 0) then
		write(*,*) "Required parameters not found in halofit."
		return
	endif

	
	!Run halofit
	status = load_matter_power(block,PK)
	if (status .ne. 0) then
		write(*,*) "Could not load matter power"
		status=3
		return
	endif

	!Load requested k range, or use same range as linear P(k)
	!if no kmin and kmax specified.
	kmin = settings%kmin
	kmax = settings%kmax
	if (kmin < 0.d0) then
		kmin = exp(PK%log_kh(1))
	endif
	if (kmax < 0.d0) then
		kmax = exp(PK%log_kh(PK%num_k))
	endif
	nk = settings%nk

	call NonLinear_GetNonLinRatios(PK)
	
	!Copy the power spectrum so we can non-linearify it
	PK_NL%num_k = PK%num_k
	PK_NL%num_z = PK%num_z
	call allocate_matterpower(PK_NL)
	PK_NL%log_kh = PK%log_kh
	PK_NL%redshifts = PK%redshifts

	!Make non-linear splines
	PK_NL%matpower = PK%matpower + 2*log(pk%nonlin_ratio)  !Since logarithmic
	call MatterPowerdata_getsplines(PK_NL)

	!Get the k values for the output - log spaced
	log_kmin = log(kmin)
	log_kmax = log(kmax)
	allocate(k(nk))
	do ik=1,nk
		k(ik) =  exp( log_kmin + (log_kmax-log_kmin)/(nk-1)*(ik-1)  )
	enddo

	!And do the interpolation
	allocate(p(nk,PK%num_z))
	do iz=1,PK%num_z
		do ik=1,nk
				p(ik,iz) = MatterPowerData_k(PK_NL,  k(ik), iz) 
				!This uses "dodgy linear interpolation" at high k.
		enddo
	enddo

	!Finally we save results
	status = datablock_put_double_grid(block, matter_power_nl_section, &
		"K_H", k, "Z", PK_NL%redshifts, "P_K", p)

	!And add another check for any errors in halofit
	if (.not. all(pk%nonlin_ratio .ge. 0)) then
		write(*,'(A)') "Halofit error.  Probably extreme parameters"
		status = 1
	endif


	deallocate(k)
	deallocate(p)
	call deallocate_matterpower(PK)
	call deallocate_matterpower(PK_NL)

end function


function cleanup(config) result(status)
	use halofit_interface_tools
	use cosmosis_modules
	type(c_ptr), value :: config
	type(halofit_settings), pointer :: settings	
	integer(cosmosis_status) :: status

	!Free memory allocated in the setup function
	call c_f_pointer(config, settings)
	deallocate(settings)

	status = 0

end function cleanup