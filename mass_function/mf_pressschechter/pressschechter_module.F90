
function setup(options) result(result)
	USE cosmosis_modules
	USE interface_tools
	implicit none
	integer(cosmosis_block), value :: options
	integer(cosmosis_status) :: status
	type(ini_settings), pointer :: settings
	type(c_ptr) :: result
	allocate(settings)
	settings%feedback = 0
        settings%redshift_zero = 0

	status = 0
	status = status + datablock_get_int_default(options, option_section, "feedback", 0, settings%feedback)
        status = status + datablock_get_int_default(options, option_section,"redshift_zero", 0, settings%redshift_zero)
	if (status .ne. 0) then 
		write(*,*) "Failed setup of PressSchechter mf!", status
		stop
	endif
	result = c_loc(settings)
end function setup

function execute(block, config) result(status)
	use compute_mf_pressschechter
	use interface_tools
	use cosmosis_modules
	implicit none
	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	type(c_ptr), value :: config
	type(ini_settings), pointer :: settings	
	type(pk_settings) :: PK
	type(massfunction) :: MassF
	real(dl), dimension(:,:), allocatable ::  dm, dr
	real(dl) :: deltam
	real(dl), allocatable, dimension(:) :: k,m
	integer iz,n

	status = 0
	call c_f_pointer(config, settings)

	

!	! mass bin sizes  in compute_mf_pressschechter.f90
	deltam = 0.01
        !In sigma.f90 lnR1=-5.684 ! 0.0034Mpc/h, 1.8e4  solar mass lnR2=4.     ! 54.9Mpc/h, 7.5e16 solar mass
        if(settings%feedback ==1) print*,"mass function (lnR1, lnR2)", lnR1, lnR2
        n =int((lnR2-lnR1)/deltam)


	
	!load in the matter power spectrum
	status = load_matter_power(block,PK)
	if (status .ne. 0) then
		write(*,*) "Could not load matter power"
		status=3
		return
	endif

        allocate(dr(n+1,PK%num_z))
        allocate(dm(n+1,PK%num_z))
        allocate(m(n+1))
        allocate(k(n+1))

        !save the z=0 output only
        if(settings%redshift_zero == 1) then
                call compute_massfunction(PK%kh,PK%matpower(:,1),MassF,n+1)
                status = datablock_put_double_array_1d(block,mass_function_section, "dndlnRh",MassF%dn_dlnRh)
                status = datablock_put_double_array_1d(block,mass_function_section, "dndlnMh",MassF%dn_dlnMh)
                status = datablock_put_double_array_1d(block,mass_function_section, "M_h",MassF%M_h)
                status = datablock_put_double_array_1d(block, mass_function_section, "R_h",MassF%R_h)
                call deallocate_mf(MassF)
        end if

        ! dr and dm = mass functions at other redshifts
        if(settings%redshift_zero == 0) then
                do iz=1,PK%num_z
                        call compute_massfunction(PK%kh,PK%matpower(:,iz),MassF,n+1)
                        k = MassF%R_h
                        m = MassF%M_h
                        dr(:,iz) = MassF%dn_dlnRh
                        dm(:,iz) = MassF%dn_dlnMh
                        call deallocate_mf(MassF)
                end do
                status =datablock_put_double_grid(block,mass_function_section,"R_H",k,"z",PK%redshifts,"dndlnRh",dr)
                status = datablock_put_double_array_1d(block, mass_function_section,"M_H",m)
                status = datablock_put_double_array_2d(block,mass_function_section,"dndlnMh", dm)
        end if



        deallocate(dr)
        deallocate(dm)
        deallocate(m)
        deallocate(k)
	call deallocate_matterpower(PK)

end function


function cleanup(config) result(status)
	use interface_tools
	use cosmosis_modules
	type(c_ptr), value :: config
	type(ini_settings), pointer :: settings	
	integer(cosmosis_status) :: status

	!Free memory allocated in the setup function
	call c_f_pointer(config, settings)
	deallocate(settings)

	status = 0

end function cleanup

