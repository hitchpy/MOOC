
program main_bvp

    use bvp_solvers, only: solve_bvp_direct
    use problem, only: u_true
    implicit none
    integer :: n, i
    real(kind=8) :: u_left, u_right, error_max, dx
    real(kind=8), dimension(:), allocatable :: x, u, ustar

    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n
    print *, "n = ",n
    allocate(x(0:n+1), u(0:n+1), ustar(0:n+1))

    dx = 1.d0 / (n+1)
    do i=0,n+1
        x(i) = i*dx
        ustar(i) = u_true(x(i))
        enddo

    u_left = 20.d0
    u_right = 60.d0

    ! check for consistency with "true solution":
    if (abs(ustar(0) - u_left) .gt. 1.d-14) then
        print *, "u_left doesn't match true solution", ustar(0)
        stop
        endif
    if (abs(ustar(n+1) - u_right) .gt. 1.d-14) then
        print *, "u_right doesn't match true solution", ustar(n+1)
        stop
        endif


    call solve_bvp_direct(x, u_left, u_right, u)


    error_max = 0.d0
    do i=1,n
        error_max = max(error_max, abs(u(i) - ustar(i)))
        enddo

    print *, "error_max = ",error_max


end program main_bvp
