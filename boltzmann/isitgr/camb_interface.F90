module camb_interface_tools
	use camb
	use cosmosis_modules
	implicit none

	character(*), parameter :: modified_gravity_section = post_friedmann_parameters_section

	integer :: standard_lmax = 1200
	integer, parameter :: CAMB_MODE_ALL = 1
	integer, parameter :: CAMB_MODE_CMB = 2
	integer, parameter :: CAMB_MODE_BG  = 3
	integer, parameter :: CAMB_MODE_THERMAL  = 4

	real(8) :: linear_zmin=0.0, linear_zmax=4.0
	integer :: linear_nz = 401

	integer :: k_eta_max_scalar = 2400
	logical :: do_lensing, do_nonlinear, do_tensors
	real(dl) :: cmb_output_scale = 7.4311e12

	real(dl), parameter :: default_yhe = 0.24
	real(dl), parameter :: default_cs2de = 1.0
	real(dl), parameter :: default_r = 0.0
	real(dl), parameter :: default_nrun = 0.0
	real(dl), parameter :: default_w = -1.0
	real(dl), parameter :: default_wa = 0.0
	real(dl), parameter :: default_pivot_scalar = 0.05
	integer,  parameter :: default_massive_nu = 0

	logical :: TGR_scale_dep= .false.  !JD for TGR 
	integer :: TGR_Rfunc = 0     !JD for TGR


	contains


	function camb_comoving_sound_horizon() result(rsdrag)
		use ModelParams
		use Precision
		use ThermoData, only : z_drag
		implicit none
		real(dl) ::  adrag, atol, rsdrag
		real(dl), external :: rombint
		integer error

		adrag = 1.0d0/(1.0d0+z_drag)
		atol = 1e-6
		rsdrag = rombint(dsound_da,1d-8,adrag,atol)
	end function camb_comoving_sound_horizon


	function camb_shift_parameter(params) result(shift_parameter)
		type(cambparams) :: params
		real(dl) :: omega_m, ombh2, omdmh2, zstar, shift_parameter
		real(dl), parameter :: c_km_per_s = 299792.458

		omega_m = params%omegac + params%omegab + params%omegan

         ombh2 = CP%omegab*(CP%h0/100.0d0)**2
         omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2

    !From Hu & Sugiyama (via modules.f90)
		zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
			(0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
			(omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

		shift_parameter = sqrt(omega_m) * params%H0 / c_km_per_s * &
		&   (1+zstar)*AngularDiameterDistance(zstar)

	end function

	function camb_initial_setup(block, mode, fixed_mode) result(status)
		integer default_lmax
		integer(c_size_t) :: block
		integer status
		character(64) :: mode_name=""
		integer :: mode
		integer, optional :: fixed_mode
		integer::  use_tabulated_w_int
		default_lmax = standard_lmax
		status=0
		! There are currently three camb modes - "background", "cmb", and "all"
		! This code may get called with a fixed mode, or with not in which case
		! we read from file
		! First in the fixed mode case, we just use that as the output mode
		if (present(fixed_mode)) then
			mode=fixed_mode
		else
			!Otherwise read from ini file

			status = datablock_get_string(block, option_section, "mode", mode_name)
			if (trim(mode_name) == "background") then
				mode=CAMB_MODE_BG
			else if (trim(mode_name) == "cmb") then
				mode=CAMB_MODE_CMB
			else if (trim(mode_name) == "all") then
				mode=CAMB_MODE_ALL
			else if (trim(mode_name) == "thermal") then
				mode=CAMB_MODE_THERMAL
			else
				write(*,*) "You need to specify a mode to use the camb module you chose."
				write(*,*) "In the camb section of your ini file, please specify one of:"
				write(*,*) "mode=background  ; For background quantities like D_A(z) only"
				write(*,*) "mode=cmb         ; For background + cmb power spectra"
				write(*,*) "mode=all         ; For background + cmb + linear matter power spectra"
				write(*,*) "mode=thermal     ; For background + thermal history params"
				write(*,*) ""
				write(*,*) "We found error status: ", status
				write(*,*) "And mode=", mode_name
				write(*,*) "Quitting now."
				stop 1
			endif
		endif

		status = 0

		!We do not use the CMB lmax if only using the background mode
		if (mode .ne. CAMB_MODE_BG) then
			status = status + datablock_get_int_default(block, option_section, "lmax", default_lmax, standard_lmax)
			status = status + datablock_get_int_default(block, option_section, "k_eta_max_scalar", 2*standard_lmax, k_eta_max_scalar)
		endif

		!We can always set an optional feedback level,
		!which defaults to zero (silent)
		status = status + datablock_get_int_default(block, option_section, "feedback", 0, FeedbackLevel)
		status = status + datablock_get_logical_default(block, option_section, "do_tensors", .false., do_tensors)

		if (mode == CAMB_MODE_ALL) then
			status = status + datablock_get_double_default(block, option_section,"zmin", linear_zmin, linear_zmin)
			status = status + datablock_get_double_default(block, option_section,"zmax", linear_zmax, linear_zmax)
			status = status + datablock_get_int_default(block, option_section,"nz", linear_nz, linear_nz)
		endif

		status = status + datablock_get_logical_default(block, option_section, "do_nonlinear", .false. , do_nonlinear)
		status = status + datablock_get_logical_default(block, option_section, "do_lensing", .false. , do_lensing)

		!Error check
		if (status .ne. 0) then
			write(*,*) "Problem setting some options for camb. Status code =  ", status
			return
		endif


 		if (do_lensing) then
 			status = status + datablock_get_string(block, option_section, "high_ell_template", highL_unlensed_cl_template)
 			if ((status .ne. 0 ) .or. trim(highL_unlensed_cl_template)=="") then
 				status = 1
 				write(*,*) "If you set do_lensing=1 then you also need to set"
 				write(*,*) "the parameter high_ell_template to the complete path"
 				write(*,*) "to the file HighLExtrapTemplate_lenspotentialCls.dat"
 				write(*,*) "which comes with CAMB - i.e. in your ini file camb section, put:"
 				write(*,*) "high_ell_template = /path/to/cosmosis/src/standard-library/boltzmann/camb/HighLExtrapTemplate_lenspotentialCls.dat"
 			elseif (.not. FileExists(trim(highL_unlensed_cl_template))) then
 				status = 2
 				write(*,*) "You set the parameter high_ell_template in the ini file to the value:"
 				write(*,*) trim(highL_unlensed_cl_template)
 				write(*,*) "But I could not find a file there.  You need to include the full"
 				write(*,*) "path to that file as that parameter, i.e.:"
 				write(*,*) "high_ell_template = /path/to/cosmosis/src/standard-library/boltzmann/camb/HighLExtrapTemplate_lenspotentialCls.dat"
 			endif
 		endif
 		
 		status = status + camb_initial_isitgr_params(block)
 		if (status .ne. 0) return
		!If noisy, report relevant params
		if (FeedbackLevel .gt. 0) then
			write(*,*) "camb mode  = ", mode
			if (mode .ne. CAMB_MODE_BG) write(*,*) "camb cmb_lmax = ", standard_lmax
			write(*,*) "camb FeedbackLevel = ", FeedbackLevel
			if (status .ne. 0) write(*,*) "Setup status: ", status
		endif
	end function camb_initial_setup

	function camb_interface_set_params(block, params, background_only) result(status)
		integer (c_int) :: status
		integer (c_size_t) :: block
		logical, optional :: background_only
		logical :: perturbations
		type(CambParams) :: params
		real(8), dimension(:), allocatable :: w_array, a_array
		character(*), parameter :: cosmo = cosmological_parameters_section
		perturbations = .true.
		if (present(background_only)) perturbations = .not. background_only

	
		call CAMB_SetDefParams(params)
		status = 0

		status = camb_set_isitgr_params(block)

        status = status + datablock_get_double(block, cosmo, "omega_b", params%omegab)
        status = status + datablock_get_double(block, cosmo, "omega_c", params%omegac)
        status = status + datablock_get_double(block, cosmo, "omega_lambda", params%omegav)
        status = status + datablock_get_double(block, cosmo, "omega_nu", params%omegan)
		status = status + datablock_get_double(block, cosmo, "omega_k", params%omegak)
		status = status + datablock_get_double(block, cosmo, "hubble", params%H0)
		
		if (perturbations) then
			status = status + datablock_get_double(block, cosmo, "n_s",     params%initpower%an(1))
            status = status + datablock_get_double_default(block, cosmo, "k_s", default_pivot_scalar, params%initpower%k_0_scalar)
			status = status + datablock_get_double(block, cosmo, "A_s",     params%initpower%ScalarPowerAmp(1))
			status = status + datablock_get_double(block, cosmo, "tau", params%Reion%optical_depth)
			status = status + datablock_get_double_default(block, cosmo, "R_T", default_r, params%initpower%rat(1))

			status = status + datablock_get_double_default(block, cosmo, "n_run", default_nrun, params%initpower%n_run(1))
			if (params%initpower%rat(1) .ne. 0) then
				status = status + datablock_get_double(block, cosmo, "n_T", params%initpower%ant(1))
			endif
		endif

		!Neutrinos

		status = status + datablock_get_double_default(block, cosmo, "cs2_de", default_cs2de, cs2_lam)
		status = status + datablock_get_double_default(block, cosmo, "yhe", default_yhe, params%yhe)

		if (params%omegan .ne. 0) then
			status = status + datablock_get_double_default(block, cosmo, "massless_nu", params%Num_Nu_massless, params%Num_Nu_massless)
			status = status + datablock_get_int_default(block, cosmo, "massive_nu", default_massive_nu, params%Num_Nu_massive)

			if (params%Num_Nu_massive == 1) then
				params%Nu_mass_eigenstates = 1
				params%nu_mass_degeneracies(1) = 1
				params%Nu_mass_fractions(1) = 1.0
			elseif (params%Num_Nu_massive == 0) then
				write(*,*) 'You need massive_nu>0 to have any omega_nu!=0'
				status=1
				return
			else
				stop "Sorry - we have not coded up the neutrino scenario your parameters implied"
			endif
		endif


		params%wantTransfer = .true.
		params%transfer%kmax = 50.0
		params%wantTensors = (params%initpower%rat(1) .ne. 0.0) .or. do_tensors

        params%Max_l=standard_lmax
        params%Max_eta_k=2*standard_lmax

        params%DoLensing = do_lensing

		!Set nonlinear behaviour
		if (do_nonlinear) then
			params%NonLinear=2
		else
			params%NonLinear=0
		endif

	
	
		!Some extras and modifications 
		params%want_zdrag = .true.
		params%want_zstar = .true.
		params%reion%use_optical_depth = .true.
		params%reion%delta_redshift = 0.5

		use_spline_template=params%DoLensing
		params%AccurateReionization = .true.
        params%transfer%num_redshifts = 1
        params%Transfer%redshifts = 0

	end function


	function camb_initial_isitgr_params(block) result (status)
		integer(cosmosis_status) ::  status
		integer(cosmosis_block) :: block
		logical use_r
		
		status = 0

		status = status + datablock_get_logical(block, &
			option_section, "scale_dependent", TGR_scale_dep)

		status = status + datablock_get_logical(block, &
			option_section, "use_r_function", use_r)

		if (status .ne. 0 ) then
			write(*,*)
			write(*,*) "Please set to T/F the variables scale_dependent and use_r_function"
			write(*,*) "In the isitgr ini file section"
			write(*,*) 
			return
		endif

        if(use_r)then
            TGR_Rfunc = 1
        else 
            TGR_Rfunc = 0
        end if

	end function

	function camb_set_isitgr_params(block) result(status)
		integer(cosmosis_status) ::  status
		integer(cosmosis_block) :: block
		real(dl) :: d_0, d_inf, r_0, r_inf, q_0, q_inf, v_0, v_inf
		
		status = 0 
		R_func = (TGR_Rfunc==1)

		status = status + datablock_get_double(block, &
			modified_gravity_section, "d_0", d_0)

		status = status + datablock_get_double(block, &
			modified_gravity_section, "d_inf", d_inf)

		status = status + datablock_get_double(block, &
			modified_gravity_section, "q_0", q_0)

		status = status + datablock_get_double(block, &
			modified_gravity_section, "q_inf", q_inf)

		status = status + datablock_get_double(block, &
			modified_gravity_section, "s", TGR_scales(2))


        if(TGR_scale_dep) then
			status = status + datablock_get_double(block, &
				modified_gravity_section, "k_c", TGR_scales(1))
		else
			TGR_scales(1) = 1e30
		endif


        r_0 = 2.d0*d_0/q_0-1.d0
        r_inf = 2.d0*d_inf/q_inf-1.d0

        v_0 = 2.d0*d_0 - q_0
        v_inf = 2.d0*d_inf - q_inf
		
		! Save these derived parameters back
		status = status + datablock_put_double(block, &
			modified_gravity_section, "v_0", v_0)

		status = status + datablock_put_double(block, &
			modified_gravity_section, "v_inf", v_inf)

		status = status + datablock_put_double(block, &
			modified_gravity_section, "r_0", r_0)

		status = status + datablock_put_double(block, &
			modified_gravity_section, "r_inf", r_inf)

		TGR_Q(1) = q_0
		TGR_Q(2) = q_inf

		if(TGR_Rfunc==1)then
		    TGR_DR(1) = r_0
		    TGR_DR(2) = r_inf
		else    
		    TGR_DR(1) = d_0
		    TGR_DR(2) = d_inf
		end if

		t_dep = (TGR_scales(2)>0.0001)

		if (status .ne. 0) then
			write(*,*) "isitgr values not set correctly"
		endif


    end function 

	
	function camb_interface_setup_zrange(params) result(status)
		integer(cosmosis_status) :: status
		type(CambParams) :: params
		real(8) :: zmin, zmax, dz
		integer nz, i

		zmin = linear_zmin
		zmax = linear_zmax
		nz = linear_nz

		dz=(zmax-zmin)/(nz-1.0)
		params%transfer%num_redshifts =  nz
        !params%Transfer%PK_num_redshifts = nz

		if (nz .gt. max_transfer_redshifts) then
			write(*,*) "Requested too many redshifts for CAMB to handle: ", nz, " = (", zmax, " - ", zmin, ") / ", dz, " + 1"
			status = 1
		endif
		
        do i=1,params%transfer%num_redshifts
			params%transfer%redshifts(nz-i+1)  = zmin + dz*(i-1)
    	enddo


    	!call Transfer_SortAndIndexRedshifts(params%transfer)
		status = 0
	end function



	function camb_interface_save_cls(block) result(status)
	
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
	
		integer, parameter :: input_set = 1
		real  :: cls(2:standard_lmax,1:4)
		real(8)  :: cls_double(2:standard_lmax,1:4), cls_phi(2:standard_lmax)
		integer  :: ell(2:standard_lmax), l
		logical, parameter :: switch_polarization_convention = .false.	
	
		status = 0
		call CAMB_GetCls(cls, standard_lmax, input_set, switch_polarization_convention)
		cls_double(:,1:4) = cls * 7.4311e12  !cmb output scale
	    do l=2,standard_lmax
			ell(l) = l
		enddo

		if (do_lensing) then
		    do l=2,standard_lmax
				cls_phi(l) = Cl_scalar(l, input_set,  C_phi) * (l+1.0) / ((l*1.0)**3 * twopi)
			enddo
		endif
	
		status = status + datablock_put_int_array_1d(block, cmb_cl_section, "ELL", ell)
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "TT", cls_double(:,1))
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "EE", cls_double(:,2))
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "BB", cls_double(:,3))
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "TE", cls_double(:,4))
		if (do_lensing) then
			status = status + datablock_put_double_array_1d(block, cmb_cl_section, "PP", cls_phi)
		endif
	
		if (status .ne. 0) then
			write(*,*) "Failed to save cmb!."
			return
		endif
	end function

	function camb_interface_save_sigma8(block) result(status)
		!Save sigma8 at z=0 to the cosmological parameters section of the file
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
		real(8) :: sigma8
		real(8), parameter :: radius8 = 8.0_8
		integer nz

		!Ask camb for sigma8
		status = 0
		sigma8=0.0
		call Transfer_Get_sigma8(MT,radius8)
		
		!It gives us the array sigma8(z).
		!We want the entry for z=0
		nz = CP%Transfer%num_redshifts
		sigma8 = MT%sigma_8(nz,1)

		!Save sigma8
		status = status + datablock_put_double(block, cosmological_parameters_section, "SIGMA_8", sigma8)
		return
	end function
	
	function camb_interface_save_transfer(block) result(status)
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
		Type(MatterPowerData) :: PK
		integer nz, nk, iz, ik
		external dtauda
		real(8) dtauda
		real(8), allocatable, dimension(:) :: k, z
		real(8), allocatable, dimension(:,:) :: P, T
		real(8), allocatable, dimension(:,:) :: ModifiedGravity_D, ModifiedGravity_D_dot
		real(8), allocatable, dimension(:,:) :: ModifiedGravity_Q, ModifiedGravity_Q_dot
		real(8) a, adotoa, TGR_D_T, TGR_D_T_dot, TGR_Q_T, TGR_Q_T_dot

		call Transfer_GetMatterPowerData(MT, PK, 1)

		nz = CP%Transfer%num_redshifts
		nk = MT%num_q_trans

		allocate(k(nk))
		allocate(z(nz))
		allocate(P(nk,nz))
		allocate(ModifiedGravity_D(nk,nz))
		allocate(ModifiedGravity_Q(nk,nz))
		allocate(ModifiedGravity_D_dot(nk,nz))
		allocate(ModifiedGravity_Q_dot(nk,nz))
		allocate(T(nk,nz))

		do ik=1,nk
			k(ik) = MT%TransferData(Transfer_kh,ik,1)
		enddo

		do iz=1,nz
			z(iz) = CP%Transfer%Redshifts(nz-iz+1)
		enddo

		do ik=1,nk
			do iz=1,nz
				P(ik,iz) = MatterPowerData_k(PK, k(ik), nz-iz+1)
				T(ik,iz) = MT%TransferData(Transfer_cdm,ik,nz-iz+1)
				!Modifications for IsItGR
				a = 1.0/(1+z(iz))
				adotoa = 1/(a*dtauda(a))
				call get_tgr_quantities(k(ik), a, adotoa, TGR_D_T, TGR_D_T_dot, TGR_Q_T, TGR_Q_T_dot)
				ModifiedGravity_D(ik,iz) = TGR_D_T
				ModifiedGravity_Q(ik,iz) = TGR_Q_T
				ModifiedGravity_D_dot(ik,iz) = TGR_D_T_dot
				ModifiedGravity_Q_dot(ik,iz) = TGR_Q_T_dot
			enddo
		enddo

		status = datablock_put_double_grid(block, matter_power_lin_section, &
        	"k_h", k, "z", z, "P_k", P)

		if (status .ne. 0) then
			write(*,*) "Failed to save matter power in CAMB."
		endif


		status = datablock_put_double_grid(block, linear_cdm_transfer_section, &
        	"k_h", k, "z", z, "delta_cdm", T)

		if (status .ne. 0) then
			write(*,*) "Failed to transfer functions in CAMB."
		endif


		status = datablock_put_double_grids(block, modified_gravity_section, &
        	"k_h", k, "z", z, "D", ModifiedGravity_D, "Q", ModifiedGravity_Q, &
        	"D_dot", ModifiedGravity_D_dot, "Q_dot", ModifiedGravity_Q_dot)

		if (status .ne. 0) then
			write(*,*) "Failed to save modified gravity functions in CAMB."
		endif


		deallocate(k, z, P, T)
	end function

	
	function camb_interface_save_da(params, block, save_density, save_thermal) result(status)
		integer (cosmosis_block) :: block
		type(CambParams) :: params
		integer (c_int) :: status
		logical, optional :: save_density, save_thermal
		logical :: density, thermal
		real(8), dimension(:), allocatable :: distance, z, rho
		character(*), parameter :: dist = distances_section
		integer nz, i
		

		! Rho as given by the code is actually 8 * pi * G * rho / c**2 , and it is measured in (Mpc)**-2
		! There really isn't a sensible set of units to do this in, so let's just use kg/m**3
		! c**2 / (8 pi G) = 5.35895884e24 kg/m
		! 1 Mpc = 3.08568025e24 m
		real(8), parameter :: mpc_in_m = 3.08568025e22
		real(8), parameter :: c2_8piG_kgm = 5.35895884e25
		real(8), parameter :: rho_units = c2_8piG_kgm / (mpc_in_m**2)
		real(8) :: shift, rs_zdrag

		density = .true.
		if (present(save_density)) density = save_density
		thermal = .true.
		if (present(save_thermal)) thermal = save_thermal

		status = 0


		nz = params%transfer%num_redshifts
		allocate(distance(nz))
		allocate(z(nz))
		if (density) allocate(rho(nz))

		do i=1,nz
			z(i) = params%transfer%redshifts(i)
			distance(i) = AngularDiameterDistance(z(i))
			!if (density) rho(i) = MT%TransferData(Transfer_rho_tot,1,i) * rho_units
		enddo
		

		shift = camb_shift_parameter(params)
		status = status + datablock_put_double(block, dist, "CMBSHIFT", shift)


!		if (thermal) then
!			status = status + datablock_put_double(block, dist, &
!				"AGE", ThermoDerivedParams( derived_Age ))
!			status = status + datablock_put_metadata(block, dist, "AGE", "unit", "Gyr")
!
!			status = status + datablock_put_double(block, dist, &
!				"RS_ZDRAG", ThermoDerivedParams( derived_rdrag ))
!
!			!There is an 
!			status = status + datablock_put_double(block, dist, &
!				"THETASTAR", ThermoDerivedParams( derived_thetastar ))
!			status = status + datablock_put_metadata(block, dist, "THETASTAR", "unit", "100 radian")
!
!			status = status + datablock_put_double(block, dist, &
!				"ZDRAG", ThermoDerivedParams( derived_zdrag ))
!
!			status = status + datablock_put_double(block, dist, &
!				"ZSTAR", ThermoDerivedParams( derived_zstar ))
!
!			status = status + datablock_put_double(block, dist, &
!				"CHISTAR", ComovingRadialDistance(ThermoDerivedParams( derived_zstar )))
!			status = status + datablock_put_metadata(block, dist, "CHISTAR", "unit", "Gyr")
!		else
			status = status + datablock_put_double(block, dist, &
				"AGE", DeltaTime(0.0_dl,1.0_dl))
			status = status + datablock_put_metadata(block, dist, "AGE", "unit", "Gyr")
!		endif


		status = status + datablock_put_double_array_1d(block, dist, "Z", z)
		status = status + datablock_put_double_array_1d(block, dist, "D_A", distance)
		status = status + datablock_put_metadata(block, dist, "D_A", "unit", "Mpc")

		distance = distance * (1+z) !Convert to D_M
		status = status + datablock_put_double_array_1d(block, dist, "D_M", distance)
		status = status + datablock_put_metadata(block, dist, "D_M", "unit", "Mpc")

		distance = distance * (1+z) !Convert to D_L
		status = status + datablock_put_double_array_1d(block, dist, "D_L", distance)
		status = status + datablock_put_metadata(block, dist, "D_L", "unit", "Mpc")

		distance = 5*log10(distance)+25 !Convert to distance modulus
		! The distance is already the dimensionful one, so we do not
		! multiply be c/H0
		status = status + datablock_put_double_array_1d(block, dist, "MU", distance)

		! Save H(z)
		!do i=1,nz
		!	distance(i) = HofZ(z(i))
		!enddo
		!status = status + datablock_put_double_array_1d(block, dist, "H", distance)
		!status = status + datablock_put_metadata(block, dist, "H", "unit", "Mpc/c")

!		if (density) then
!			status = status + datablock_put_double_array_1d(block, dist, "RHO", rho)
!			status = status + datablock_put_metadata(block, dist, "RHO", "unit", "KG/M^3")
!		endif


		status = status + datablock_put_int(block, dist, "NZ", nz)

		!And finally save a
		z = 1.0/(1+z)
		status = status + datablock_put_double_array_1d(block, dist, "A", z)


		if (status .ne. 0) then
			write(*,*) "Failed to write redshift-distance column data in block section."
		endif
		
		deallocate(distance)
		deallocate(z)
!		if (density) deallocate(rho)
		
	end function
	
end module camb_interface_tools



