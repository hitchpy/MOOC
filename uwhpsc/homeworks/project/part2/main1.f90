
program main_heat1

    use heat_solvers, only: solve_heat_explicit
    use problem, only: k, u_true
    implicit none
    integer :: n, i, nsteps
    real(kind=8) :: error_max, dx, pi, tfinal, t0
    real(kind=8), dimension(:), allocatable :: x, u, ustar, u0

    pi = acos(-1.d0)

    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n
    read(21,*) k
    read(21,*) tfinal
    read(21,*) nsteps
    print '("n = ",i6)', n
    print '("k = ",i6)', k
    print '("tfinal = ",f7.4)', tfinal
    print '("nsteps = ",i6)', nsteps

    allocate(x(0:n+1), u(0:n+1), ustar(0:n+1), u0(0:n+1))

    t0 = 0.d0
    dx = pi / (n+1)
    do i=0,n+1
        x(i) = i*dx
        enddo


    do i=0,n+1
        u0(i) = sin(k*x(i))
        ustar(i) = u_true(x(i), tfinal)
        enddo


    call solve_heat_explicit(x, u0, t0, tfinal, nsteps, u)


    error_max = 0.d0
    do i=0,n+1
        ustar(i) = u_true(x(i), tfinal)
        error_max = max(error_max, abs(u(i) - ustar(i)))
        enddo

    print '("error_max = ", e13.6)', error_max

    open(unit=22, file='solution.txt', status='unknown')
    do i=0,n+1
        write(22,220) x(i), u(i), ustar(i)
220     format(3e22.14)
        enddo

    close(22)

end program main_heat1
