! Liklehood from massively simplified WMAP data
! Just looks at the CMB shift parameter R
! Reported by WMAP9 as 
! R = 1.728 +/- 0.016

function execute(block) result(status)
	use cosmosis_modules
	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	real(8) :: R, like

	!Fixed parameter values for R_mu and R_sigma
	real(8), parameter ::  R_mu = 1.728
	real(8), parameter ::  R_sigma = 0.016
	
	!Error checking
	status = 0

	!Extract the shift parameter and calculate the likelihood
	R = 0.0
	status = datablock_get_double(block, distances_section, "cmbshift", R)
	like = -0.5*(R - R_mu)**2/R_sigma**2

	!Save the likelihood and close the file
	status = status + datablock_put_double(block, likelihoods_section, "SHIFT_LIKE", like)

	!status is the return value - any problems will be passed along

end function