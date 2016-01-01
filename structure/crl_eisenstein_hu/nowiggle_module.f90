
function setup(options) result(result)
  USE cosmosis_modules
  USE interface_tools
  implicit none
  integer(cosmosis_block), value :: options
  integer(cosmosis_status) :: status
  type(ini_settings), pointer :: settings
  type(c_ptr) :: result
  allocate(settings)

  status = 0
  status = status + datablock_get_double_default(options, option_section, "zmin", 0.0D-04, settings%zmin)
  status = status + datablock_get_double_default(options, option_section, "zmax", 3.0D+00, settings%zmax)
  status = status + datablock_get_int_default(options, option_section, "nz_steps", 800, settings%nz_steps)
  status = status + datablock_get_double_default(options, option_section, "kmin", 1.0D-05, settings%kmin)
  status = status + datablock_get_double_default(options, option_section, "kmax", 10D+00, settings%kmax)
  status = status + datablock_get_int_default(options, option_section, "nk_steps", 800, settings%nk_steps)

  settings%dz=(settings%zmax-settings%zmin)/(settings%nz_steps-1.0)
  settings%dk=(log(settings%kmax)-log(settings%kmin))/(settings%nk_steps-1.0)

  print*,""
  print*,"EH power spectrum options:"
  print*, "zmin",settings%zmin
  print*, "zmax",settings%zmax
  print*, "nz_steps",settings%nz_steps
  print*, "kmin",settings%kmin
  print*, "kmax",settings%kmax
  print*, "nk_steps",settings%nk_steps
  print*,""

  if (status .ne. 0) then
     write(*,*) "Failed setup of CRL POWER SPECTRUM NO WIGGLE!", status
     stop
  endif
  result = c_loc(settings)
end function setup



! ==========================================

function execute(block, config) result(status)
  use compute_pk_nowiggle
  use interface_tools
  use cosmosis_modules
  implicit none
  integer(cosmosis_block), value :: block
  integer(cosmosis_status) :: status
  type(c_ptr), value :: config
  type(ini_settings), pointer :: settings
  type(pk_settings) :: PK
  real(8) :: omega_baryon,omega_matter,w_de,omega_de,h0,n_s,n_run,A_s
  real(8) :: k, z,D_0,yval, ypval, yppval
  real(8), allocatable, dimension(:,:) :: P
  real(8), allocatable, dimension(:) :: dz,zbins,dz_interpolated


  integer iz,n,i,n_growth

  status = 0
  call c_f_pointer(config, settings)


  !  LOAD COMSMOLOGICAL PARAMETERS

  status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_B", omega_baryon)
  status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_M", omega_matter)
  status = status + datablock_get_double_default(block, cosmological_parameters_section, "W",-1.0d0, w_de)
  status = status + datablock_get_double(block, cosmological_parameters_section, "h0", h0)
  status = status + datablock_get_double(block, cosmological_parameters_section, "n_s", n_s)
  status = status + datablock_get_double_default(block, cosmological_parameters_section, "n_run", 0.0d0, n_run)
  status = status + datablock_get_double(block, cosmological_parameters_section, "A_s", A_s)
  n_growth= datablock_get_array_length(block, GROWTH_PARAMETERS_SECTION, "d_z")
  status = status+ datablock_get_double_array_1d(block, GROWTH_PARAMETERS_SECTION, "d_z", dz,n_growth);
  status = status+  datablock_get_double_array_1d(block, GROWTH_PARAMETERS_SECTION, "z", zbins,n_growth);

  if (status .ne. 0) then
    write(*,*) "Error in crl_eistenstein_hu"
    return
  endif

  !    INTERPOLATE GROWTH
  allocate(dz_interpolated(n_growth))

  call spline_cubic_set ( n_growth, zbins , dz, 2, 0.0, 2, 0.0, dz_interpolated )

  if ( (abs(zbins(1) - settings%zmin) .gt. 0.1d0) .or. abs(zbins(n_growth)- settings%zmax).gt. 0.1d0 ) then

     print*, "======================="
     print*,  " the chosen bounds of the growth module does not cover the entire redshift range requested for the power spectrum"
     print*, "zmin zmax from growth",zbins(1),zbins(n_growth)
     print*, "zmin zmax in Pk_nowiggle",settings%zmin,settings%zmax

     print*, "======================="
  end if



  ! create an array for P(k,z)

  ! first define the dimensions
  PK%num_k = settings%nk_steps
  PK%num_z = settings%nz_steps

  ! create empty containers, defined in interface.tools.f90

  call  allocate_matterpower(PK)

  z=settings%zmin
  k=log(settings%kmin)

  !  fill k and z arrays

  do i = 1, PK%num_z, 1
     PK%redshifts(i)= z
     z =z+settings%dz
  end do

  do i = 1, PK%num_k, 1

     PK%kh(i)= exp(k)
     k =k+settings%dk

  end do

  ! finally fill with P(k,z) with EH formula
  ! cycle over k


  do i = 1, PK%num_k, 1
     k=PK%kh(i)
     call   compute_pknowiggle(k,A_s,n_s,h0,omega_matter,omega_baryon,PK%matpower(i,1))
  end do
  ! cycle over z just applying the growth

  D_0=dz_growth(PK%redshifts(1),zbins,dz,dz_interpolated)
  PK%matpower(:,1)=PK%matpower(:,1)*(D_0)**2
  do i = 2, PK%num_z, 1
     z=PK%redshifts(i)
     PK%matpower(:,i)=PK%matpower(:,1)*(dz_growth(z,zbins,dz,dz_interpolated)/D_0)**2
  end do


  ! save in datablock

  status = datablock_put_double_grid(block, "matter_power_no_bao", &
       "k_h", PK%kh, "z", PK%redshifts, "P_k", PK%matpower)


  !! deallocate everything here.

  call deallocate_matterpower(PK)

  if (allocated(dz_interpolated)) deallocate(dz_interpolated)

end function execute






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

