module mead_settings_mod
	type mead_settings
		logical :: noisy
		real(8) :: kmin, kmax
		integer :: nk

		real(8) :: numin, numax

		real(8) :: zmin, zmax
		integer :: nz

		logical :: feedback

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

	status = status + datablock_get(options, option_section, "zmin", settings%zmin)
	status = status + datablock_get(options, option_section, "zmax", settings%zmax)
	status = status + datablock_get(options, option_section, "nz", settings%nz)


	status = status + datablock_get(options, option_section, "kmin", settings%kmin)
	status = status + datablock_get(options, option_section, "kmax", settings%kmax)
	status = status + datablock_get(options, option_section, "nk", settings%nk)

	status = status + datablock_get_double_default(options, option_section, "numin", 0.1D0, settings%numin)
	status = status + datablock_get_double_default(options, option_section, "numax", 5.0D0, settings%numax)

	status = status + datablock_get_logical_default(options, option_section, "feedback", .false., settings%feedback)

	if (status .ne. 0) then
		write(*,*) "One or more parameters not found for hmcode"
		stop
	endif

	WRITE(*,*) 'z min:', settings%zmin
	WRITE(*,*) 'z max:', settings%zmax
	WRITE(*,*) 'number of z:', settings%nz
	WRITE(*,*)

	WRITE(*,*) 'k min:', settings%kmin
	WRITE(*,*) 'k max:', settings%kmax
	WRITE(*,*) 'number of k:', settings%nk
	WRITE(*,*)


	result = c_loc(settings)

end function setup


function execute(block,config) result(status)
	use mead_settings_mod
	use cosmosis_modules
	use mhm

	implicit none

	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	type(c_ptr), value :: config
	type(mead_settings), pointer :: settings	
	integer, parameter :: LINEAR_SPACING = 0
	integer, parameter :: LOG_SPACING = 1
	character(*), parameter :: cosmo = cosmological_parameters_section
	character(*), parameter :: halo = halo_model_parameters_section
	character(*), parameter :: linear_power = matter_power_lin_section
	character(*), parameter :: nl_power = matter_power_nl_section

	real(4) :: p1h, p2h,pfull, plin, z
	integer :: i,j, z_index
	REAL, ALLOCATABLE :: k(:),  pmod(:), ztab(:)
	TYPE(cosmology) :: cosi
	TYPE(tables) :: lut
	!CosmoSIS supplies double precision - need to convert
	real(8) :: om_m, om_v, om_b, h, w, sig8, n_s, om_nu
	real(8), ALLOCATABLE :: k_in(:), z_in(:), p_in(:,:)
	real(8), ALLOCATABLE :: k_out(:), z_out(:), p_out(:,:)
	real(8) :: Halo_as, halo_eta0

	status = 0
	call c_f_pointer(config, settings)

	feedback = settings%feedback

	!Fill in the cosmology parameters. We need to convert from CosmoSIS 8-byte reals
	!to HMcode 4-byte reals, hence the extra bit
	status = status + datablock_get(block, cosmo, "omega_m", om_m)
	status = status + datablock_get(block, cosmo, "omega_lambda", om_v)
	status = status + datablock_get(block, cosmo, "omega_b", om_b)
        status = status + datablock_get_double_default(block, cosmo, "omega_nu", 0.0D0, om_nu)
	status = status + datablock_get(block, cosmo, "h0", h)
	status = status + datablock_get(block, cosmo, "sigma_8", sig8)
	status = status + datablock_get(block, cosmo, "n_s", n_s)
	status = status + datablock_get_double_default(block, cosmo, "w", -1.0D0, w)


	status = status + datablock_get_double_default(block, halo, "A", 3.13D0, halo_as)
	status = status + datablock_get_double_default(block, halo, "eta_0", 0.603D0, halo_eta0)

	if (status .ne. 0 ) then
		write(*,*) "Error reading parameters for Mead code"
		return
	endif

    cosi%om_m=om_m-om_nu !The halo modelling should include only cold matter components (i.e. DM and baryons)
    cosi%om_v=om_v
    cosi%om_b=om_b
    cosi%h=h
    cosi%w=w
    cosi%sig8=sig8
    cosi%n=n_s

    cosi%eta_0 = halo_eta0
    cosi%As = halo_as

    !And get the cosmo power spectrum, again as double precision
    !Also the P is 2D as we get z also
	status = status + datablock_get_double_grid(block, linear_power, &
        "k_h", k_in, "z", z_in, "p_k", p_in)

	if (status .ne. 0 ) then
		write(*,*) "Error reading P(k,z) for Mead code"
		return
	endif

	!Copy in k
	allocate(cosi%ktab(size(k_in)))
	cosi%ktab = k_in

	!Find the index of z where z==0
	if (z_in(1)==0.0) then
		z_index=1
	elseif (z_in(size(z_in))==0.0) then
		z_index=size(z_in)
	else
		write(*,*) "P(k,z=0) not found - please calculate"
		status = 1
		return
	endif
	!Copy in P(k) from the right part of P(k,z)
	allocate(cosi%pktab(size(k_in)))
    cosi%pktab = p_in(:, z_index) * (cosi%ktab**3.)/(2.*(pi**2.))
    cosi%itk = 5


	!Set the output ranges in k and z
	CALL fill_table(real(settings%kmin),real(settings%kmax),k,settings%nk,LOG_SPACING)
	CALL fill_table(real(settings%zmin),real(settings%zmax),ztab,settings%nz,LINEAR_SPACING)

	!Fill table for output power
	ALLOCATE(p_out(settings%nk,settings%nz))


	!Loop over redshifts
	DO j=1,settings%nz

		!Sets the redshift
		z=ztab(j)

		!Initiliasation for the halomodel calcualtion
		!Also normalises power spectrum (via sigma_8)
		!and fills sigma(R) tables
		CALL halomod_init(z,real(settings%numin),real(settings%numax),lut,cosi)

		!Loop over k values
		DO i=1,SIZE(k)
			plin=p_lin(k(i),cosi)        
			CALL halomod(k(i),z,p1h,p2h,pfull,plin,lut,cosi)
			!This outputs k^3 P(k).  We convert back.
			p_out(i,j)=pfull / (k(i)**3.0) * (2.*(pi**2.))
		END DO

		IF(j==1) THEN
			if (settings%feedback) WRITE(*,fmt='(A5,A7)') 'i', 'z'
			if (settings%feedback) WRITE(*,fmt='(A13)') '   ============'
		END IF
		 if (settings%feedback) WRITE(*,fmt='(I5,F8.3)') j, ztab(j)
	END DO

	!convert to double precision
	allocate(k_out(settings%nk))
	allocate(z_out(settings%nz))
	k_out = k
	z_out = ztab
	!Convert k to k/h to match other modules
	!Output results to cosmosis
	status = datablock_put_double_grid(block,nl_power, "k_h", k_out, "z", z_out, "p_k", p_out)

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
    IF(ALLOCATED(cosi%rtab)) DEALLOCATE(cosi%rtab)
    IF(ALLOCATED(cosi%sigtab)) DEALLOCATE(cosi%sigtab)   
    IF(ALLOCATED(cosi%ktab)) DEALLOCATE(cosi%ktab)
    IF(ALLOCATED(cosi%tktab)) DEALLOCATE(cosi%tktab)
    IF(ALLOCATED(cosi%pktab)) DEALLOCATE(cosi%pktab)

end function execute
