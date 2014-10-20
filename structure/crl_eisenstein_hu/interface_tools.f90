MODULE interface_tools
  use cosmosis_modules
  IMPLICIT none




  integer, parameter :: dl=8
  type ini_settings
     real(dl):: zmin,zmax,kmin,kmax,dz,dk
     integer ::nk_steps ,nz_steps
  end type ini_settings

  type pk_settings
     real(dl), dimension(:), allocatable :: redshifts
     real(dl), dimension(:), allocatable :: kh
     real(dl), dimension(:,:), allocatable :: matpower
     real(dl), dimension(:,:), allocatable :: ddmat
     integer num_z, num_k
  end type pk_settings



contains

  function load_matter_power(block, PK) result(status)
    use cosmosis_modules
    integer(cosmosis_block) :: block
    integer(cosmosis_status) :: status
    type(pk_settings) :: PK
    real(dl), allocatable, dimension(:) :: k, z
    real(dl), allocatable, dimension(:,:) :: P

    !Get the data columns from the fits data
    status = 0

    !Load k, z, P
    status = datablock_get_double_grid(block, matter_power_lin_section, &
         "K_H", k, "Z", z, "P_K", P)

    if (status .ne. 0) then
       write(*,*) "Could not find K_H, Z, or P_K in block"
       return
    endif

    !Fill in data structure
    PK%num_k = size(k)
    PK%num_z = size(z)
    call allocate_matterpower(PK)
    PK%kh = k
    PK%redshifts = z
    PK%matpower = P


    !Clean up
    deallocate(k, z, P)



  end function load_matter_power



    DOUBLE PRECISION FUNCTION dz_growth(z,xa,ya,ya2) ! growth function fitting
    IMPLICIT none
    DOUBLE PRECISION, intent(IN) :: z
    DOUBLE PRECISION, allocatable, dimension(:) , intent(IN) :: xa,ya,ya2
    DOUBLE PRECISION :: a,b,hh,x, ypval, yppval
    INTEGER :: jlo,dim
    dim=0
    dim=size(ya, dim=1)
    x=z+0.00000000000000000001d0
    call spline_cubic_val( dim, xa, ya, ya2, z, dz_growth, ypval, yppval)
    return
  END FUNCTION dz_growth



  subroutine allocate_matterpower(PK)
    type(pk_settings) :: PK
    allocate(PK%redshifts(PK%num_z))
    allocate(PK%kh(PK%num_k))
    allocate(PK%matpower(PK%num_k,PK%num_z))
    allocate(PK%ddmat(PK%num_k,PK%num_z))
  end subroutine allocate_matterpower

  subroutine deallocate_matterpower(PK)
    type(pk_settings) :: PK
    deallocate(PK%redshifts)
    deallocate(PK%kh)
    deallocate(PK%matpower)
    deallocate(PK%ddmat)
  end subroutine deallocate_matterpower


END MODULE interface_tools
