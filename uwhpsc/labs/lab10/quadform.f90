! Sample solution to the Lab 10 problem of computing a quadratic form

program quadform

    use omp_lib
    implicit none
    integer, parameter :: n=10
    integer :: nthreads, thread_num

    real (kind=8), dimension(n,n) :: A 
    real (kind=8), dimension(n) :: x, Ax
    real (kind=8) :: xAx

    integer :: i,j

    nthreads = 1  ! if not using OpenMP
    !$ call omp_set_num_threads(6)


    ! initialize matrix A and vector x:
    A = 0.d0
    !$omp parallel do 
    do i=1,n
        A(i,i) = 1.d0
        x(i) = i
        enddo

    ! multiply A*x to get Ax:
    !$omp parallel do private(j)
    do i=1,n
        Ax(i) = 0.d0
        do j=1,n
            Ax(i) = Ax(i) + A(i,j)*x(j)
            enddo
        enddo

    xAx = 0.d0
    !$omp parallel do reduction(+:xAx)
    do i=1,n
        xAx = xAx + x(i) * Ax(i)
        enddo

    print 600, xAx
600 format("xAx = ", e22.15)

end program quadform
