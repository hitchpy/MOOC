
program hilbert

    implicit none
    real(kind=8), dimension(:),allocatable :: x,b,work
    real(kind=8), dimension(:,:),allocatable :: a
    real(kind=8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i, info, lda, ldb, nrhs, n, j
    integer, dimension(:), allocatable :: ipiv
    integer, allocatable, dimension(:) :: iwork
    character, dimension(1) :: norm
    logical :: debug

    debug = .false.

    print *, "Input n... "
    read *, n

    allocate(a(n,n))
    allocate(b(n))
    allocate(x(n))
    allocate(ipiv(n))

    nrhs = 1 ! number of right hand sides in b
    lda = n  ! leading dimension of a
    ldb = n  ! leading dimension of b

    do j=1,n
        do i=1,n
            a(i,j) = 1.d0 / (i+j-1.d0)
            enddo
        enddo

    x = 1.d0
    b = matmul(a,x)

    
    if (debug) then
        print *, "A = "
        do i=1,n
            write(6,601) (a(i,j), j=1,n)
601         format(20e13.6)
            enddo
        endif

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    do i=1,n
        print *, b(i)
        enddo

    if (debug) then
        print *, "A factors = "
        do i=1,n
            write(6,601) (a(i,j), j=1,n)
            enddo
        endif
  

    ! compute 1-norm of error
    errnorm = 0.d0
    xnorm = 0.d0
    do i=1,n
        errnorm = errnorm + abs(x(i)-b(i))
        xnorm = xnorm + abs(x(i))
        enddo

    ! relative error in 1-norm:
    errnorm = errnorm / xnorm

    ! compute 1-norm needed for condition number

    anorm = 0.d0
    do j=1,n
        colsum = 0.d0
        do i=1,n
            colsum = colsum + abs(a(i,j))
            enddo
        anorm = max(anorm, colsum)
        enddo

    ! compute condition number of matrix:
    ! note: uses A returned from dgesv with L,U factors:

    allocate(work(4*n))
    allocate(iwork(n))
    norm = '1'  ! use 1-norm
    call dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)

    if (info /= 0) then
        print *, "*** Error in dgecon: info = ",info
        endif

    print *, 'anorm = ',anorm
    print *, 'rcond = ',rcond
    print 201, n, 1.d0/rcond, errnorm
201 format("For n = ",i4," the approx. condition number is ",e12.3,/, &
           " and the relative error in 1-norm is ",e10.3)    


    deallocate(a,b,ipiv)
    deallocate(work,iwork)

end program hilbert
