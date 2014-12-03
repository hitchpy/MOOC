
program main

    use random_util2, only: init_random_seed

    implicit none
    integer :: seed
    real (kind=8) :: x(3)

    print *, "input seed"
    read *, seed

    call init_random_seed(seed)  

    call random_number(x)
    print *, "Three random numbers: ",x
    call random_number(x)
    print *, "Three more random numbers: ",x


end program main

