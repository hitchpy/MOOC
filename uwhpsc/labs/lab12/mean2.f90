
! Sample solution to the problem posed in Lab 12.
! Version with coarse grain parallelism

program main

    use omp_lib
    use random_util, only: init_random_seed

    implicit none
    integer :: seed, nthreads
    integer, parameter :: ntests = 7
    real (kind=8), allocatable :: x(:)
    real (kind=8) :: sum, count1, count4
    real (kind=8), dimension(ntests)  :: mean,q1,q4
    integer :: i,n,k,kstart,kend,thread_num,points_per_thread

    ! Specify number of threads to use:
    nthreads = 1       ! need this value in serial mode
    !$ nthreads = 2
    !$ call omp_set_num_threads(nthreads)
    !$ print "('Using OpenMP with ',i3,' threads')", nthreads

    print *, "input seed"
    read *, seed

    call init_random_seed(seed)  

    print *, "      n          mean        quartile 1    quartile 4"

    points_per_thread = (ntests + nthreads - 1) / nthreads
    print *, "points_per_thread = ",points_per_thread

    !$omp parallel private(k,x,count1,count4,sum,kstart,kend,thread_num,n)

    thread_num = 0     ! needed in serial mode
    !$ thread_num = omp_get_thread_num()    ! unique for each thread

    ! Determine start and end index for the set of points to be
    ! handled by this thread:
    kstart = thread_num * points_per_thread + 1
    kend = min((thread_num+1) * points_per_thread, ntests)

    !$omp critical
    print 201, thread_num, kstart, kend
    !$omp end critical
201 format("Thread ",i2," will take k = ",i6," through k = ",i6)

    do k=kstart,kend
        n = 10**(k+1)
        allocate(x(n))
        call random_number(x)
        sum = 0.d0
        count1 = 0
        count4 = 0
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
        enddo
    !$omp end parallel


    do k=1,ntests
        n = 10**(k+1)
        print 601, n,mean(k),q1(k),q4(k)
 601    format(i10,3f15.8)
        enddo
    

end program main

