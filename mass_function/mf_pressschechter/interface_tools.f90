MODULE interface_tools
  use cosmosis_modules
  IMPLICIT none
        integer, parameter :: dl=8
        type ini_settings
                integer :: feedback
                integer :: redshift_zero
        end type
        
        type massfunction 
                real(dl), dimension(:), allocatable :: M_h
                real(dl), dimension(:), allocatable :: R_h
                real(dl), dimension(:), allocatable :: dn_dlnRh
                real(dl), dimension(:), allocatable :: dn_dlnMh
                integer num_r, num_z
        end type
                

        type pk_settings
                real(dl), dimension(:), allocatable :: redshifts
                real(dl), dimension(:), allocatable :: kh
                real(dl), dimension(:,:), allocatable :: matpower
                real(dl), dimension(:,:), allocatable :: ddmat
                integer num_z, num_k
        end type

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

        !Allocate memory
        !call MatterPowerdata_getsplines(PK)

        !Clean up
        deallocate(k, z, P)

end function


                subroutine allocate_matterpower(PK)
                        type(pk_settings) :: PK
                        allocate(PK%redshifts(PK%num_z))
                        allocate(PK%kh(PK%num_k))
                        allocate(PK%matpower(PK%num_k,PK%num_z))
                        allocate(PK%ddmat(PK%num_k,PK%num_z))
                end subroutine

                subroutine deallocate_matterpower(PK)
                        type(pk_settings) :: PK
                        deallocate(PK%redshifts)
                        deallocate(PK%kh)
                        deallocate(PK%matpower)
                        deallocate(PK%ddmat)
                end subroutine

                subroutine allocate_mf(MassF)
                        type(massfunction) :: MassF
                        allocate(MassF%M_h(MassF%num_r))
                        allocate(MassF%R_h(MassF%num_r))
                        allocate(MassF%dn_dlnRh(MassF%num_r))
                        allocate(MassF%dn_dlnMh(MassF%num_r))
                end subroutine

                subroutine deallocate_mf(MassF)
                        type(massfunction) :: MassF
                        deallocate(MassF%M_h)
                        deallocate(MassF%R_h)
                        deallocate(MassF%dn_dlnRh)
                        deallocate(MassF%dn_dlnMh)
                end subroutine

END MODULE interface_tools
