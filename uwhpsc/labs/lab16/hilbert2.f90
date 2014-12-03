program hilbert2

    implicit none
    integer :: n
    real(kind=8) :: cond

    open(21, file='cond.txt',status='unknown')

    do n=2,20
        cond = hilbert_condition(n)
        !print *, "cond = ",cond
        write(21, 210) n,cond
 210    format(i4,e16.6)
        enddo

    print *, "Created cond.txt"


contains

    real(kind=8) function hilbert_condition(n)

    implicit none
    integer, intent(in) :: n
    real(kind=8), dimension(:),allocatable :: x,b,work
    real(kind=8), dimension(:,:),allocatable :: a
    real(kind=8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i, info, lda, ldb, nrhs, j
    integer, dimension(:), allocatable :: ipiv
    integer, allocatable, dimension(:) :: iwork
    character, dimension(1) :: norm

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

    ! do not need to solve a system, but do need to do factorization:
    call dgetrf(n, n, a, lda, ipiv, info)

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

    hilbert_condition = 1.d0 / rcond

    deallocate(a,b,ipiv)
    deallocate(work,iwork)

    end function hilbert_condition

end program hilbert2
