module ModelParams
	implicit None
	!This is just the camb params used in halofit
	integer, parameter :: dl = 8
    real(dl), parameter :: pi = 3.141592654

	type CambParams
		real(dl) :: omegac
		real(dl) :: omegab
		real(dl) :: omegan
		real(dl) :: omegav		
	end type CambParams

	type(CambParams) :: CP
end module

module LambdaGeneral
    use  ModelParams
    implicit none

    real(dl)  :: w_lam = -1_dl !p/rho for the dark energy (assumed constant)
    ! w_lam is now w0
    !comoving sound speed. Always exactly 1 for quintessence
    !(otherwise assumed constant, though this is almost certainly unrealistic)
    real(dl) :: cs2_lam = 1_dl
    !cs2_lam now is ce^2

    logical :: use_tabulated_w = .false.
    real(dl) :: wa_ppf = 0._dl
    real(dl) :: c_Gamma_ppf = 0.4_dl
    integer, parameter :: nwmax = 5000, nde = 2000
    integer :: nw_ppf
    real(dl) w_ppf(nwmax), a_ppf(nwmax), ddw_ppf(nwmax)
    real(dl) rde(nde),ade(nde),ddrde(nde)
    real(dl), parameter :: amin = 1.d-9
    logical :: is_cosmological_constant
    private nde,ddw_ppf,rde,ade,ddrde,amin
end module LambdaGeneral