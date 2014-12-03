
! Sample solution to the problem posed in Lab 12.

program main

    use omp_lib
    use random_util, only: init_random_seed

    implicit none
    integer :: seed, nthreads
    integer, parameter :: ntests = 7
    real (kind=8), allocatable :: x(:)
    real (kind=8) :: sum, count1, count4
    real (kind=8), dimension(ntests)  :: mean,q1,q4
    integer :: i,n,k

    !$ nthreads = 2
    !$ call omp_set_num_threads(nthreads)
    !$ print *, "Number of threads: ", nthreads

    print *, "input seed"
    read *, seed

    call init_random_seed(seed)  

    print *, "      n          mean        quartile 1    quartile 4"
    do k=1,ntests
        n = 10**(k+1)
        allocate(x(n))
        call random_number(x)
        sum = 0.d0
        count1 = 0
        count4 = 0
        !$omp parallel do reduction(+:sum,count1,count4)
        do i=1,n
            sum = sum + x(i)
            if (x(i) < 0.25) then
                count1 = count1 + 1
                endif
            if (x(i) > 0.75) then
                count4 = count4 + 1
                endif
            enddo
        mean(k) = sum / float(n)
        q1(k) = count1 / float(n)
        q4(k) = count4 / float(n)
        deallocate(x)

        print 601, n,mean(k),q1(k),q4(k)
 601    format(i10,3f15.8)
        enddo
    

end program main

