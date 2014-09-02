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
	use ModelParams
	implicit None
	real(dl) :: w_lam
end module