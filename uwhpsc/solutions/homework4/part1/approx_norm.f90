

module approx_norm

    implicit none
    integer :: nsamples, seed
    save

contains

    subroutine approx_norm1(a, anorm)

    use omp_lib
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: a
    real(kind=8), intent(out) :: anorm

    
    real(kind=8), allocatable, dimension(:) :: x,ax,allx
    integer :: i,j,k,n,k1,k2
    integer :: kstart,kend,thread_num,nthreads,samples_per_thread
    real(kind=8) :: xnorm, axnorm, ratio, max_ratio, max_ratio_thread

    n = size(a,1) 
    allocate(x(n),ax(n))

    ! generate all the random numbers we'll ever need:
    allocate(allx(n*nsamples))
    call random_number(allx)



    max_ratio = 0.d0

    !$omp parallel private(i,j,max_ratio_thread,x,ax,xnorm,axnorm, &
    !$omp                  kstart,kend,thread_num,samples_per_thread) 

    thread_num = 0     ! needed in serial mode
    !$ thread_num = omp_get_thread_num()    ! unique for each thread

    nthreads = 1  ! needed in serial mode
    !$ nthreads = omp_get_num_threads()
    samples_per_thread = (nsamples + nthreads - 1) / nthreads

    ! Determine start and end index for the set of samples to be 
    ! handled by this thread:
    kstart = thread_num * samples_per_thread + 1
    kend = min((thread_num+1) * samples_per_thread, nsamples)

    !$omp critical
    print 201, thread_num, kstart, kend
    !$omp end critical
201 format("Thread ",i2," will take k = ",i6," through k = ",i6)

    max_ratio_thread = 0.d0

    do k=kstart,kend

        ! choose a random vector x by taking n elements from allx :
        k1 = (k-1)*n+1
        k2 = k*n
        x = allx(k1:k2)

        ! compute matrix-vector product ax:
        do i=1,n    
            ax(i) = 0.d0
            do j=1,n
                ax(i) = ax(i) + a(i,j)*x(j)
                enddo
            enddo

        ! compute 1-norm of x and ax
        xnorm = 0.d0
        axnorm = 0.d0
        do i=1,n
            xnorm = xnorm + abs(x(i))
            axnorm = axnorm + abs(ax(i))
            enddo

        ratio = axnorm/xnorm
        max_ratio_thread = max(max_ratio_thread, ratio)
        enddo
    

    ! update global max_ratio with value from each thread:
    !$omp critical
      max_ratio = max(max_ratio, max_ratio_thread)
    !$omp end critical

    !$omp end parallel 

    anorm = max_ratio

end subroutine approx_norm1

end module approx_norm
