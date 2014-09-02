module Transfer
	use ModelParams
	implicit none


    Type MatterPowerData
        !everything is a function of k/h
        integer   ::  num_k, num_z
        real(dl), dimension(:), pointer :: log_kh => NULL(), redshifts => NULL()
        !matpower is log(P_k)
        real(dl), dimension(:,:), allocatable :: matpower, ddmat 
        !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
        !function of k and redshift NonLinearScaling(k_index,z_index)
        real(dl), dimension(:,:), pointer :: nonlin_ratio => NULL()
    end Type MatterPowerData
	contains


    subroutine MatterPowerdata_getsplines(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

    do i = 1,PK_Data%num_z
        call spline(PK_data%log_kh,PK_data%matpower(1,i),PK_data%num_k,&
        cllo,clhi,PK_data%ddmat(1,i))
    end do

	end subroutine MatterPowerdata_getsplines


    function MatterPowerData_k(PK,  kh, itf) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: kh
    real(dl) :: logk
    integer llo,lhi
    real(dl) outpower, dp
    real(dl) ho,a0,b0
    integer, save :: i_last = 1

    logk = log(kh)
    if (logk < PK%log_kh(1)) then
        dp = (PK%matpower(2,itf) -  PK%matpower(1,itf)) / &
        ( PK%log_kh(2)-PK%log_kh(1) )
        outpower = PK%matpower(1,itf) + dp*(logk - PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter

        dp = (PK%matpower(PK%num_k,itf) -  PK%matpower(PK%num_k-1,itf)) / &
        ( PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1) )
        outpower = PK%matpower(PK%num_k,itf) + dp*(logk - PK%log_kh(PK%num_k))
    else
        llo=min(i_last,PK%num_k)
        do while (PK%log_kh(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_kh(llo+1)< logk)
            llo=llo+1
        end do
        i_last =llo
        lhi=llo+1
        ho=PK%log_kh(lhi)-PK%log_kh(llo)
        a0=(PK%log_kh(lhi)-logk)/ho
        b0=1-a0

        outpower = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
        ((a0**3-a0)* PK%ddmat(llo,itf) &
        +(b0**3-b0)*PK%ddmat(lhi,itf))*ho**2/6
    end if

    outpower = exp(max(-30._dl,outpower))

    end function MatterPowerData_k

	SUBROUTINE spline(x,y,n,yp1,ypn,y2)
	implicit none
	INTEGER, intent(in) :: n
	real(dl), intent(in) :: x(n), y(n), yp1, ypn
	real(dl), intent(out) :: y2(n)
	INTEGER i,k
	real(dl) p,qn,sig,un
	real(dl), dimension(:), allocatable :: u


	Allocate(u(1:n))
	if (yp1.gt..99d30) then
		y2(1)=0._dl
		u(1)=0._dl
	else
		y2(1)=-0.5d0
		u(1)=(3._dl/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	endif

	do i=2,n-1
		sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
		p=sig*y2(i-1)+2._dl 

		y2(i)=(sig-1._dl)/p

		u(i)=(6._dl*((y(i+1)-y(i))/(x(i+ &
		1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
		u(i-1))/p
	end do
	if (ypn.gt..99d30) then
		qn=0._dl
		un=0._dl
	else
		qn=0.5d0
		un=(3._dl/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._dl)
	do k=n-1,1,-1
		y2(k)=y2(k)*y2(k+1)+u(k)
	end do

	Deallocate(u)

	!  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
	END SUBROUTINE spline
	
	subroutine allocate_matterpower(PK)
		type(MatterPowerData) :: PK
		allocate(PK%redshifts(PK%num_z))
		allocate(PK%log_kh(PK%num_k))
		allocate(PK%matpower(PK%num_k,PK%num_z))
		allocate(PK%ddmat(PK%num_k,PK%num_z))
	    allocate(PK%nonlin_ratio(PK%num_k,PK%num_z))
	end subroutine
	
	subroutine deallocate_matterpower(PK)
		type(MatterPowerData) :: PK
		deallocate(PK%redshifts)
		deallocate(PK%log_kh)
		deallocate(PK%matpower)
		deallocate(PK%ddmat)
    	deallocate(PK%nonlin_ratio)
	end subroutine

end module Transfer