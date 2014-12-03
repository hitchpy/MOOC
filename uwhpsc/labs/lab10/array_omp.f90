
program array1

    implicit none
    integer, parameter :: m = 5, n=6
    integer :: nthreads, thread_num

    real (kind=8), dimension(m,n) :: A 
    real (kind=8), dimension(m) :: b
    real (kind=8), dimension(n) :: x 

    real (kind=8), dimension(m,n) :: At
    real (kind=8), dimension(m) :: bt 
    real (kind=8), dimension(n) :: xt 
    integer :: i,j

    nthreads = 1  ! if not using OpenMP
    !$ call omp_set_num_threads(6)
    !$ nthreads = omp_get_num_threads()
    print *, "nthreads = ",nthreads 
 
    ! initialize matrix A and vector x:
    do i=1,m
        $omp parallel do private(i)
        do j=1,n
            A(i,j) = i+j
            enddo
        enddo

    ! serial version for At:
    do i=1,m
        do j=1,n
            At(i,j)=i+j
            enddo
        enddo

    xt = 1.
    x  = 1.

    ! multiply A*x to get b:
    !$omp parallel do private(j)
    do i=1,m
        b(i) = 0.
        do j=1,n
            b(i) = b(i) + A(i,j)*x(j)
            enddo
        enddo

    bt=matmul(At,xt)

    print *, "b and bt should be equal "
    print *, "b="
    print "(d16.6)", b
    print *, "bt="
    print "(d16.6)",bt

end program array1
