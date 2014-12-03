
module heat_solvers

contains

subroutine solve_heat_explicit(x, u0, t0, tfinal, nsteps, u)
    implicit none
    real(kind=8), intent(in) :: x(0:), u0(0:), t0, tfinal
    integer, intent(in) :: nsteps
    real(kind=8), intent(out) :: u(0:)
    integer :: n, i, nstep
    real(kind=8) :: dx, dt, dtdx2, um1, up1
    real(kind=8), allocatable :: uxx(:)

    n = size(x) - 2
    dx = x(2) - x(1)
    dt = (tfinal - t0) / float(nsteps)

    print 21, n, dx, dt
21  format("Solving heat equation with n = ",i6, " dx = ",e10.3," dt = ",e10.3)

    dtdx2 = dt / dx**2
    print 22, dtdx2
22  format("  dt/dx**2 = ",f7.3)
	if (dtdx2 > 0.5d0) then
		print *, "*** Warning: the explicit method is unstable"
		endif

    allocate(uxx(1:n))

    u = u0

    do nstep=1,nsteps

        do i=1,n
            uxx(i) = (u(i-1) - 2.d0*u(i) + u(i+1)) / dx**2
            enddo

        do i=1,n
            u(i) = u(i) + dt * uxx(i)
            enddo

        enddo

end subroutine solve_heat_explicit


end module heat_solvers
