! This sample solution to Homework 3 has only one method,
! to use as a starting point for Homework 4, part 1.

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
    real(kind=8) :: xnorm, axnorm, ratio, max_ratio

    n = size(a,1) 
    allocate(x(n),ax(n))

    ! generate all the random numbers we'll ever need:
    allocate(allx(n*nsamples))
    call random_number(allx)

    max_ratio = 0.d0

    !$omp parallel do private(i,j,x,ax,xnorm,axnorm,ratio) &
    !$omp reduction(max:max_ratio)
    do k=1,nsamples

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
        max_ratio = max(max_ratio, ratio)
        enddo

    anorm = max_ratio

    end subroutine approx_norm1

end module approx_norm
