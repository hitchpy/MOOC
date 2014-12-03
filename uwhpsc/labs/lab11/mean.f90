
! Sample solution to the problem posed in Lab 11.

program main

    use omp_lib
    use random_util2, only: init_random_seed

    implicit none
    integer :: seed
    real (kind=8), allocatable :: x(:)
    real (kind=8) :: mean
    integer :: i,n

    print *, "input n"
    read *, n
    allocate(x(n))
    print *, "input seed"
    read *, seed

    call init_random_seed(seed)  

    call random_number(x)
    mean = 0.d0
    !$omp parallel do reduction(+:mean)
    do i=1,n
        mean = mean + x(i)
        enddo
    mean = mean / n
    print *, "The mean is: ",mean

    print *, "The difference from the true mean 0.5 is ",abs(mean - 0.5d0)
    print *, "1/sqrt(n) = ", 1.d0/sqrt(float(n))


end program main

