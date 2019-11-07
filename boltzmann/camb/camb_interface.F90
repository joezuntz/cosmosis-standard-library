module camb_interface_tools
	use camb
	use cosmosis_modules
	implicit none

	integer :: standard_lmax = 1200
	real(dl) :: standard_kmax = 50.0
	!Not the actual values camb uses; may involve
	!extrapolation
	real(8) :: standard_kmin
	integer :: standard_nk
	integer, parameter :: CAMB_MODE_ALL = 1
	integer, parameter :: CAMB_MODE_CMB = 2
	integer, parameter :: CAMB_MODE_BG  = 3
	integer, parameter :: CAMB_MODE_THERMAL  = 4

	real(8) :: linear_zmin=0.0, linear_zmax=4.0
	integer :: linear_nz = 101
	
	!Vivian begins
	real(8) :: back_zmin=0.0, back_zmax=4.0
	integer :: back_nz = 401
	!Vivian ends

	logical :: distances_to_lss
	integer :: n_highz_distance

	integer :: k_eta_max_scalar = 2400
	logical :: do_lensing, do_nonlinear, do_tensors

	integer :: matter_power_lin_version = 1  ! Consider this a bit-field
	                                 ! indicating which spectra are to
	                                 ! be produced.
	integer, parameter :: MATTER_POWER_LIN_            =  1
	integer, parameter :: MATTER_POWER_LIN_CDM_BARYON  =  2
	integer, parameter :: MATTER_POWER_LIN_BOTH        =  3
	integer, parameter :: MATTER_POWER_LIN_WEYL        =  4
	integer, parameter :: MATTER_POWER_LIN_WEYL_MATTER =  5
	integer, parameter :: MATTER_POWER_LIN_ALL         =  8

	integer :: sterile_neutrino = 0
	real(dl) :: delta_neff = 0.0
	real(dl) :: sterile_mass_fraction = 0.0
	real(dl) :: cmb_output_scale = 7.4311e12

	real(dl), parameter :: default_yhe = 0.24
	real(dl), parameter :: default_cs2de = 1.0
	real(dl), parameter :: default_r = 0.0
	real(dl), parameter :: default_nrun = 0.0
	real(dl), parameter :: default_w = -1.0
	real(dl), parameter :: default_wa = 0.0
	real(dl), parameter :: default_pivot_scalar = 0.05
	integer,  parameter :: default_massive_nu = 0
	integer,  parameter :: default_sterile_neutrinos = 0
	real(dl),  parameter :: default_kmax = 50.0
	real(dl),  parameter :: default_kmin = -1.0 ! default is to use whatever camb this is sensible
	integer, parameter :: default_nk = -1  ! default is to use whatever camb this is sensible

	logical :: need_thermal_init

	!From Vivian Miranda, to expose accuracy settings
	real(dl) :: acc_boost = 1.0_dl
	logical :: high_acc_default = .false.

	logical :: reject_non_accelerating
	logical :: save_growth


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

		need_thermal_init = (mode==CAMB_MODE_THERMAL)

		!We do not use the CMB lmax if only using the background mode
		if (mode .ne. CAMB_MODE_BG) then
			status = status + datablock_get_int_default(block, option_section, "lmax", default_lmax, standard_lmax)
			status = status + datablock_get_int_default(block, option_section, "k_eta_max_scalar", 2*standard_lmax, k_eta_max_scalar)
		endif

		!We can always set an optional feedback level,
		!which defaults to zero (silent)
		status = status + datablock_get_int_default(block, option_section, "feedback", 0, FeedbackLevel)
		status = status + datablock_get_logical_default(block, option_section, "use_tabulated_w", .false., use_tabulated_w)
		status = status + datablock_get_logical_default(block, option_section, "do_tensors", .false., do_tensors)

		status = status + datablock_get_double_default(block, option_section,"zmin", linear_zmin, linear_zmin)
		status = status + datablock_get_double_default(block, option_section,"zmax", linear_zmax, linear_zmax)
		status = status + datablock_get_int_default(block, option_section,"nz", linear_nz, linear_nz)

		status = status + datablock_get_double_default(block, option_section,"background_zmin", back_zmin, back_zmin)
		status = status + datablock_get_double_default(block, option_section,"background_zmax", back_zmax, back_zmax)
		status = status + datablock_get_int_default(block, option_section,"background_nz", back_nz, back_nz)


		status = status + datablock_get_double_default(block, option_section,"kmax", default_kmax, standard_kmax)
		status = status + datablock_get_double_default(block, option_section,"kmin", default_kmin, standard_kmin)
		status = status + datablock_get_int_default(block, option_section,"nk", default_nk, standard_nk)

		!Either neither or both of nk, kmin must be set
		if ((standard_nk<=0 .and. standard_kmin>0) .or. (standard_nk>0 .and. standard_kmin<=0)) then
			write(*,*) "If you set either kmin or nk for camb then you must set both (and to a positive number)."
			status = 1
		endif


		status = status + datablock_get_logical_default(block, option_section, "do_nonlinear", .false. , do_nonlinear)
		status = status + datablock_get_logical_default(block, option_section, "do_lensing", .false. , do_lensing)
		status = status + datablock_get_int_default(block, option_section, "matter_power_lin_version", MATTER_POWER_LIN_, matter_power_lin_version)


		status = status + datablock_get_logical_default(block, option_section, "distances_to_lss", .false. , distances_to_lss)
		status = status + datablock_get_int_default(block, option_section, "n_highz_distance", 100, n_highz_distance)

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
 		
	    !VM BEGINS
	    status = status + datablock_get_double_default(block,option_section,"accuracy_boost",1.0_dl, acc_boost)
	    AccuracyBoost = acc_boost
	    status = status + datablock_get_logical_default(block,option_section,"high_accuracy_default",.false.,high_acc_default)
	    HighAccuracyDefault = high_acc_default

	    status = status + datablock_get_logical_default(block, option_section, "reject_non_accelerating", .false., reject_non_accelerating)
	    status = status + datablock_get_logical_default(block, option_section, "save_growth", .false., save_growth)

		!If noisy, report relevant params
		if (FeedbackLevel .gt. 0) then
			write(*,*) "camb mode  = ", mode
			if (mode .ne. CAMB_MODE_BG) write(*,*) "camb cmb_lmax = ", standard_lmax
			write(*,*) "camb FeedbackLevel = ", FeedbackLevel
			if (status .ne. 0) write(*,*) "Setup status: ", status
			write(*,*) "accuracy boost = ", AccuracyBoost
			write(*,*) "HighAccuracyDefault = ", HighAccuracyDefault
		endif
	end function camb_initial_setup

	function camb_interface_set_params(block, params, mode) result(status)
		integer (c_int) :: status
		integer (c_size_t) :: block
		integer :: mode
		logical :: perturbations
		type(CambParams) :: params
		integer :: sterile_neutrino
		real(8) :: nu_mass_1
		character(9) :: param_name
		integer :: i
		real(8), dimension(:), allocatable :: w_array, a_array
		character(*), parameter :: cosmo = cosmological_parameters_section
		perturbations = (mode .eq. CAMB_MODE_CMB) .or. (mode .eq. CAMB_MODE_ALL)

	
		call CAMB_SetDefParams(params)
		status = 0

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
			status = status + datablock_get_int_default(block, cosmo, "sterile_neutrino", default_sterile_neutrinos, sterile_neutrino)
			status = status + datablock_get_double_default(block, cosmo, "massless_nu", params%Num_Nu_massless, params%Num_Nu_massless)
			status = status + datablock_get_int_default(block, cosmo, "massive_nu", default_massive_nu, params%Num_Nu_massive)
			! We check to see 
			status = status + datablock_get_double_default(block, cosmo, "nu_mass_1", -1000.0d0, nu_mass_1)

			!  Scenario 1: sterile neutrinos
			if (sterile_neutrino > 0) then
				status = status + datablock_get_double(block, cosmo, "delta_neff", delta_neff)
				status = status + datablock_get_double(block, cosmo, "sterile_mass_fraction", sterile_mass_fraction)
				params%share_delta_neff = .false.
				params%Num_Nu_massless = 2.0307
				params%Nu_mass_eigenstates = 2
				params%Num_Nu_massive = 2
				params%nu_mass_degeneracies(1) = 1.0153
				params%nu_mass_degeneracies(2) = delta_neff
				params%nu_mass_fractions(1) = (1.0 - sterile_mass_fraction) 
				params%nu_mass_fractions(2) = sterile_mass_fraction
			! Scenario 2: single massive nu
			elseif (params%Num_Nu_massive == 1) then
				params%Nu_mass_eigenstates = 1
				params%nu_mass_numbers(1) = 1
				params%Nu_mass_fractions(1) = 1.0
				params%share_delta_neff = .true.
			! Scenario 3: individually specified masses
			elseif (nu_mass_1 .ge. -999.) then
				! Loop through the neutrinos that we have and
				! get a mass for each
				params%nu_mass_numbers(1:params%Num_Nu_massive) = 1
				params%Nu_mass_eigenstates = params%Num_Nu_massive
				do i=1,params%Num_Nu_massive
					write(param_name, '(A8,I1)') "nu_mass_", i
					status = status + datablock_get_double_default(block, cosmo, param_name, -1000.0d0, params%Nu_mass_fractions(i))
					if (params%Nu_mass_fractions(i) < -999.) then
						write(*,*) "Must set all nu masses if you set any of them"
						status = status + 1
						exit
					endif

				enddo
				params%Nu_mass_fractions(1:params%Num_Nu_massive) = &
								params%Nu_mass_fractions(1:params%Num_Nu_massive) / &
								sum(params%Nu_mass_fractions(1:params%Num_Nu_massive))
				params%share_delta_neff = .true.
			! Scenario 3: one mass eigenstate
			elseif (params%Num_Nu_massive == 3) then
				params%Nu_mass_eigenstates = 1
				params%nu_mass_numbers(1) = 3
				params%Nu_mass_fractions(1) = 1.0
				params%share_delta_neff = .true.
			elseif (params%Num_Nu_massive == 0) then
				write(*,*) 'You need massive_nu>0 to have any omega_nu!=0'
				status=1
				return
			else
				stop "Sorry - we have not coded up the neutrino scenario your parameters implied"
			endif
		endif


		call setcgammappf()

		! tabulated dark energy EoS
		if (use_tabulated_w) then
			status = status + datablock_get_double_array_1d(block, de_equation_of_state_section, "w", w_array, nw_ppf)
			status = status + datablock_get_double_array_1d(block, de_equation_of_state_section, "a", a_array, nw_ppf)
			if (status .ne. 0) then
				write(*,*) ""
				write(*,*) "CAMB TABULATED W(A) ERROR:"
				write(*,*) "Failed to read w(a) from de_equation_of_state"
				write(*,*) ""
				return
			endif
			if (nw_ppf .gt. nwmax) then
				write(*,*) ""
				write(*,*) "CAMB TABULATED W(A) ERROR:"
				write(*,*) "The size of the w(a) table was too large ", nw_ppf, nwmax
				write(*,*) ""
				status=nw_ppf
				return
			endif
			w_ppf(1:nw_ppf) = w_array(1:nw_ppf)
			a_ppf(1:nw_ppf) = dlog(a_array(1:nw_ppf))  !a is stored as log(a)
			if (a_array(1) > a_array(2)) then
				write(*,*) ""
				write(*,*) "CAMB TABULATED W(A) ERROR:"
				write(*,*) "The tabulated w(a) in camb should have the scale factor a in increasing order"
				write(*,*) ""
				status = 1
				return
			endif
			deallocate(w_array, a_array)
			call setddwa()
			call interpolrde()
		else
			status = status + datablock_get_double_default(block, cosmo, "w", -1.0D0, w_lam)
			status = status + datablock_get_double_default(block, cosmo, "wa",  0.0D0, wa_ppf)
			if (reject_non_accelerating .and. (w_lam+wa_ppf .gt. 0)) then
				write(*,*) "Unphysical w_0 + w_a = ", w_lam, " + ", wa_ppf, " = ", w_lam+wa_ppf, " > 0"
				status = 1
			endif
		endif	


		params%wantTransfer = (mode==CAMB_MODE_ALL)
		params%transfer%kmax = standard_kmax
		params%wantTensors = (params%initpower%rat(1) .ne. 0.0) .or. do_tensors

        params%Max_l=standard_lmax
        params%Max_eta_k=2*standard_lmax

        params%DoLensing = do_lensing
        params%DerivedParameters = .true.
		!Set nonlinear behaviour
		if (do_nonlinear) then
			params%NonLinear=2
		else
			params%NonLinear=0
		endif

		if (mode==CAMB_MODE_THERMAL) params%reion%Reionization = .false.
	
		!Some extras and modifications 
		params%want_zdrag = .true.
		params%want_zstar = .true.
		params%reion%use_optical_depth = .true.
		params%reion%delta_redshift = 0.5

		use_spline_template=params%DoLensing
		params%AccurateReionization = .true.
        params%Transfer%PK_num_redshifts = 1
        params%Transfer%PK_redshifts = 0

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
		params%transfer%num_redshifts = nz
        params%Transfer%PK_num_redshifts = nz

		if (nz .gt. max_transfer_redshifts) then
			write(*,*) "Requested too many redshifts for CAMB to handle: ", nz, " = (", zmax, " - ", zmin, ") / ", dz, " + 1"
			status = 1
		endif
		
        do i=1,params%transfer%num_redshifts
			params%transfer%redshifts(nz-i+1)  = zmin + dz*(i-1)
	        params%transfer%pk_redshifts(nz-i+1)  = zmin + dz*(i-1)
    	enddo

    	call Transfer_SortAndIndexRedshifts(params%transfer)
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
	


	function  camb_interface_save_transfer (block)  result(status)
   		integer (cosmosis_block)  ::  block
		integer  (cosmosis_status)  ::  status

        	status  =  0

		if  (iand (matter_power_lin_version, MATTER_POWER_LIN_)  .ne.  0)  then
			status = status + camb_interface_save_transfer__ (block, MATTER_POWER_LIN_)
		endif

		if  (iand (matter_power_lin_version, MATTER_POWER_LIN_CDM_BARYON)  .ne.  0)  then
			status  =  status  +  camb_interface_save_transfer__ (block, MATTER_POWER_LIN_CDM_BARYON)
		endif

		if  (iand (matter_power_lin_version, MATTER_POWER_LIN_WEYL)  .ne.  0)  then
			status  =  status  +  camb_interface_save_transfer__ (block, MATTER_POWER_LIN_WEYL)
		endif

      end function camb_interface_save_transfer



	function camb_interface_save_transfer__ (block, matter_power_lin_version__)  result(status)
		integer (cosmosis_block) :: block
		integer  ::  matter_power_lin_version__
		integer (cosmosis_status) :: status
		Type(MatterPowerData) :: PK
		integer nz, nk, iz, ik

		real(8), allocatable, dimension(:) :: k, z
		real(8), allocatable, dimension(:,:) :: T

		if  (matter_power_lin_version__  .eq.  MATTER_POWER_LIN_)   then
			call Transfer_GetMatterPowerData  (MT, PK, 1)
		elseif (matter_power_lin_version__  .eq.  MATTER_POWER_LIN_CDM_BARYON)  then
			call Transfer_GetMatterPowerData  (MT, PK, 1, var1=transfer_nonu, var2=transfer_nonu)
		elseif (matter_power_lin_version__  .eq.  MATTER_POWER_LIN_WEYL)  then
			call Transfer_GetMatterPowerData  (MT, PK, 1, var1=transfer_Weyl, var2=transfer_Weyl)
		endif

		nz = CP%Transfer%num_redshifts
		nk = MT%num_q_trans

		allocate(k(nk))
		allocate(z(nz))
		allocate(T(nk,nz))

		do ik=1,nk
			k(ik) = MT%TransferData(Transfer_kh,ik,1)
		enddo

		do iz=1,nz
			z(iz) = CP%Transfer%Redshifts(nz-iz+1)
		enddo

		do ik=1,nk
			do iz=1,nz
				T(ik,iz) = MT%TransferData(Transfer_cdm,ik,nz-iz+1)
			enddo
		enddo


        	status = 0

        	if (.not.  datablock_has_section (block, linear_cdm_transfer_section))  then
			status = datablock_put_double_grid(block, linear_cdm_transfer_section, &
						"k_h", k, "z", z, "delta_cdm", T)

			if (status .ne. 0) then
				write(*,*) "Failed to save transfer function in CAMB."
			endif

		endif

		deallocate(k, z, T)

		!Now save the matter power
		status = status + camb_interface_save_matter_power(block, PK, matter_power_lin_version__)

		call MatterPowerdata_Free(PK)
	end function



	function camb_interface_save_matter_power(block, PK, matter_power_lin_version__)  result(status)
		integer (cosmosis_block) :: block
		integer  ::  matter_power_lin_version__
		integer (cosmosis_status) :: status
		Type(MatterPowerData) :: PK
		integer nz, nk, iz, ik
		real(8) :: log_kmin, log_kmax
		real(8), allocatable, dimension(:) :: k, z
		real(8), allocatable, dimension(:,:) :: P
		character(50) :: datablock_section

		! z values to use.  Use sample values from camb,
		! which were set in the input
		nz = CP%Transfer%num_redshifts
		allocate(z(nz))
		do iz=1,nz
			z(iz) = CP%Transfer%Redshifts(nz-iz+1)
		enddo


		!k values.  Two cases. Either the user set in the
		!input file and we use what they asked for, or they didn't
		!and we just use whatever defaults camb chose.
		!This latter one is for backwards compatibility.
		if (standard_nk>0) then
			nk = standard_nk			
			log_kmin = log(standard_kmin)
			log_kmax = log(standard_kmax)
			allocate(k(nk))
			do ik=1,nk
				k(ik) = exp(log_kmin + (log_kmax-log_kmin)*(ik-1)/(nk-1))
			enddo
		else
			nk = MT%num_q_trans
			allocate(k(nk))
			do ik=1,nk
				k(ik) = MT%TransferData(Transfer_kh,ik,1)
			enddo
		endif


		allocate(P(nk,nz))


		do ik=1,nk
			do iz=1,nz
				P(ik,iz) = MatterPowerData_k(PK, k(ik), nz-iz+1)
			enddo
		enddo

		if  (matter_power_lin_version__  .eq.  MATTER_POWER_LIN_)  then
			datablock_section  =  matter_power_lin_section
		elseif  (matter_power_lin_version__  .eq.  MATTER_POWER_LIN_CDM_BARYON)  then
			datablock_section  =  matter_power_lin_cdm_baryon_section
		elseif  (matter_power_lin_version__  .eq.  MATTER_POWER_LIN_WEYL)  then
			datablock_section  =  "weyl_curvature_spectrum"
		endif

		status = datablock_put_double_grid(block, datablock_section, &
					"k_h", k, "z", z, "P_k", P)

		if (status .ne. 0) then
			write(*,*) "Failed to save matter power in CAMB."
		endif

		deallocate(k, z, P)

	end function


	
	function camb_interface_save_da(params, block, save_density, save_thermal) result(status)
		use CAMBmain
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
		real(8), parameter :: z_lss = 1200.0
		real(8) dlogz_high, zmax_regular
		integer nz_regular


		density = .true.
		if (present(save_density)) density = save_density
		thermal = .true.
		if (present(save_thermal)) thermal = save_thermal

		status = 0


!Vivian begins
		!nz = params%transfer%num_redshifts		
		!nz_regular = nz
		nz = back_nz
		nz_regular = nz
!Vivian ends
		!Fixed number of sample points out to last scattering surface
		if (distances_to_lss) then 
			nz = nz + n_highz_distance
!Vivian begins			
			!zmax_regular = params%transfer%redshifts(1)
			zmax_regular = back_zmax
!Vivian ends			
			dlogz_high = (log(1+z_lss) - log(1+zmax_regular))/n_highz_distance
		endif

		allocate(distance(nz))
		allocate(z(nz))
		!if (density) allocate(rho(nz))

!Vivian begins
		do i=1,nz
			if (i<=nz_regular) then
				! The low redshift regime
				!z(i) = params%transfer%redshifts(nz_regular-i+1)
				z(i) = back_zmin + (i-1)*(back_zmax-back_zmin)/(back_nz-1.0)
!Vivian ends				
			else
				! Additional points beyond the regime where we store the matter power spectrum
				z(i) = (1+zmax_regular) * exp(dlogz_high*(i-nz_regular)) - 1
			endif
			distance(i) = AngularDiameterDistance(z(i))
			!if (density) rho(i) = MT%TransferData(Transfer_rho_tot,1,i) * rho_units
		enddo


		!Need to call thermal history here
		if (thermal .and. need_thermal_init) then
			call cmbmain()
		endif

		shift = camb_shift_parameter(params)
		status = status + datablock_put_double(block, dist, "CMBSHIFT", shift)


		if (thermal) then
			status = status + datablock_put_double(block, dist, &
				"AGE", ThermoDerivedParams( derived_Age ))
			status = status + datablock_put_metadata(block, dist, "AGE", "unit", "Gyr")

			status = status + datablock_put_double(block, dist, &
				"RS_ZDRAG", ThermoDerivedParams( derived_rdrag ))
			status = status + datablock_put_metadata(block, dist, "RS_ZDRAG", "unit", "Mpc")

			!There is an 
			status = status + datablock_put_double(block, dist, &
				"THETASTAR", ThermoDerivedParams( derived_thetastar ))
			status = status + datablock_put_metadata(block, dist, "THETASTAR", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"ZDRAG", ThermoDerivedParams( derived_zdrag ))


			status = status + datablock_put_double(block, dist, &
				"K_D", ThermoDerivedParams( derived_kD ))
			status = status + datablock_put_metadata(block, dist, "K_D", "unit", "1/Mpc")

			status = status + datablock_put_double(block, dist, &
				"THETA_D", ThermoDerivedParams( derived_thetaD ))
			status = status + datablock_put_metadata(block, dist, "THETA_D", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"Z_EQUALITY", ThermoDerivedParams( derived_zEQ ))

			status = status + datablock_put_double(block, dist, &
				"K_EQUALITY", ThermoDerivedParams( derived_keq ))
			status = status + datablock_put_metadata(block, dist, "K_EQUALITY", "unit", "1/Mpc")


			status = status + datablock_put_double(block, dist, &
				"THETA_EQUALITY", ThermoDerivedParams( derived_thetaEQ ))
			status = status + datablock_put_metadata(block, dist, "THETA_EQUALITY", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"THETA_RS_EQUALITY", ThermoDerivedParams( derived_theta_rs_EQ ))
			status = status + datablock_put_metadata(block, dist, "THETA_RS_EQUALITY", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"DA_STAR", ThermoDerivedParams( derived_DAstar ))
			status = status + datablock_put_metadata(block, dist, "DA_STAR", "unit", "Gpc")

			status = status + datablock_put_double(block, dist, &
				"R_STAR", ThermoDerivedParams( derived_rstar ))
			status = status + datablock_put_metadata(block, dist, "R_STAR", "unit", "Mpc")

			status = status + datablock_put_double(block, dist, &
				"ZSTAR", ThermoDerivedParams( derived_zstar ))

			status = status + datablock_put_double(block, dist, &
				"CHISTAR", ComovingRadialDistance(ThermoDerivedParams( derived_zstar )))
			status = status + datablock_put_metadata(block, dist, "CHISTAR", "unit", "Mpc")

			status = status + datablock_put_double(block, dist, &
				"THETA_MC", CosmomcTheta())
			status = status + datablock_put_metadata(block, dist, "THETA_MC", "unit", "radian")

		else
			status = status + datablock_put_double(block, dist, &
				"AGE", DeltaPhysicalTimeGyr(0.0_dl,1.0_dl))
			status = status + datablock_put_metadata(block, dist, "AGE", "unit", "Gyr")
		endif


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
		do i=1,nz
			distance(i) = HofZ(z(i))
		enddo
		status = status + datablock_put_double_array_1d(block, dist, "H", distance)
		status = status + datablock_put_metadata(block, dist, "H", "unit", "1.0/Mpc")

		!if (density) then
		!	status = status + datablock_put_double_array_1d(block, dist, "RHO", rho)
		!	status = status + datablock_put_metadata(block, dist, "RHO", "unit", "KG/M^3")
		!endif


		status = status + datablock_put_int(block, dist, "NZ", nz)

		!And finally save a
		z = 1.0/(1+z)
		status = status + datablock_put_double_array_1d(block, dist, "A", z)


		if (status .ne. 0) then
			write(*,*) "Failed to write redshift-distance column data in block section."
		endif
		
		deallocate(distance)
		deallocate(z)
		!if (density) deallocate(rho)
		
	end function

	function camb_interface_save_growth_rate(block) result(status)
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
		real(8), allocatable, dimension(:) :: z, sigma8_z, sigma2_vdelta_8_z, fsigma8_z, f_z, d_z
		integer nz, iz
		real(8), parameter :: radius8 = 8.0_8

		status = 0

		if (.not. save_growth) return

		!Ask camb for sigmas
		call Transfer_Get_sigmas(MT,radius8)

		! collect info into arrays we're going to save
		nz = CP%Transfer%num_redshifts
		allocate(z(nz))
		allocate(sigma8_z(nz))
		allocate(sigma2_vdelta_8_z(nz))
		allocate(fsigma8_z(nz))
		allocate(f_z(nz))
		allocate(d_z(nz)) 

		do iz=1,nz
			z(iz) = CP%Transfer%Redshifts(nz-iz+1)
			sigma8_z(iz) = MT%sigma_8(nz-iz+1,1)
			sigma2_vdelta_8_z(iz) = MT%sigma2_vdelta_8(nz-iz+1,1)
			fsigma8_z(iz) = sigma2_vdelta_8_z(iz)/sigma8_z(iz)             
			f_z(iz) = fsigma8_z(iz)/sigma8_z(iz)
			d_z(iz) = sigma8_z(iz)/sigma8_z(1)
		enddo

		! f*sigma8(z) is defined in terms of power spectra as (sigma_vd_8(z))^2/sigma8(z),
		!   where sigma_vd_8(z) measures the smoothed density-velocity correlations defined
		!   anologously to sigma8, but using the velocity-density power spectrum P_vd, where v is
		!   the Newtonian-gauge peculier velocity of baryons and DM, while d is the total matter
		!   density fluctuation (following Planck cosmology papers)

		!Save to datablock
		!to do set up name for growth secction
		status = status + datablock_put_double_array_1d(block, growth_parameters_section, "SIGMA_8_Z", sigma8_z)
		status = status + datablock_put_double_array_1d(block, growth_parameters_section, "SIGMA2_VDELTA_8_Z", sigma2_vdelta_8_z)
		status = status + datablock_put_double_array_1d(block, growth_parameters_section, "FSIGMA8_Z", fsigma8_z)
		status = status + datablock_put_double_array_1d(block, growth_parameters_section, "F_Z", f_z)
		status = status + datablock_put_double_array_1d(block, growth_parameters_section, "D_Z", d_z) 
		status = status + datablock_put_double_array_1d(block, growth_parameters_section, "Z", z)

		return
	end function camb_interface_save_growth_rate

	
end module camb_interface_tools
