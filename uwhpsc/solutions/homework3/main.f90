! Timing statements have been added to this version for illustration purposes.

program main

    use omp_lib
    use random_util, only: init_random_seed
    use approx_norm, only: approx_norm1, nsamples, method

    implicit none
    real(kind=8), allocatable, dimension(:,:) :: a
    real(kind=8) :: colsum, colsum_max, anorm, t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate
    integer :: n, seed, i, j, nthreads, seed_matrix

    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n
    read(21,*) seed
    read(21,*) nsamples
    read(21,*) nthreads
    read(21,*) method

    !$ call omp_set_num_threads(nthreads)
    print *, "nthreads = ",nthreads

    call init_random_seed(seed)

    ! Create random matrix:
    allocate(a(n,n))
    call random_number(a)


    ! Compute 1-norm as max column sum:
    colsum_max = 0.d0
    !$omp parallel do private(i,colsum) reduction(max: colsum_max)
    do j=1,n
        colsum = 0.d0
        do i=1,n
            colsum = colsum + abs(a(i,j))
            enddo
        colsum_max = max(colsum_max, colsum)
        enddo
            
    !---- check the time for the call to approx_norm:

    call system_clock(tclock1)
    call cpu_time(t1)   ! start cpu timer

    ! Estimate 1-norm:
    call approx_norm1(a, anorm)

    call cpu_time(t2)   ! end cpu timer
    print 10, t2-t1
 10 format("CPU time = ",f12.8, " seconds")

    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")


    print 600, n, colsum_max
600 format("Using matrix with n = ",i5,"   True 1-norm:   ",f16.6)
    print 601, nsamples, anorm
601 format("Based on ",i6," samples, Approximate 1-norm: ",f16.6)

end program main
