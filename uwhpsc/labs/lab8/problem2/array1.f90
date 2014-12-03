! $UWHPSC/codes/fortran/array1

program array1

    ! demonstrate declaring and using arrays

    implicit none
    integer, parameter :: m = 5, n=6

    real (kind=8), dimension(n,m) :: A 
    real (kind=8), dimension(m) :: b
    real (kind=8), dimension(n) :: x 

    real (kind=8), dimension(m,n) :: At
    real (kind=8), dimension(m) :: bt 
    real (kind=8), dimension(n) :: xt 
    integer :: i,j

    ! initialize matrix A and vector x:
    do i=1,m
        do j=1,n
            A(i,j) = i+j
            At(i,j)=i+j
            enddo
        enddo

    xt = 1.
    x  = 1.

    ! multiply A*x to get b:
    do i=1,m
        b(i) = 0.
        do j=1,m
            b(i) = b(i) + A(i,j)*x(j)
            enddo
        enddo
    bt=matmul(At,xt)
    print *, "A = "
    do i=1,m
        print *, A(i,:)   ! i'th row of A
        enddo
    print "(2d16.6)", ((A(i,j), j=1,2), i=1,3)
    print *, "x = "
    print "(d16.6)", x
    print *, "b and bt should be equal "
    print *, "b="
    print "(d16.6)", b
    print *, "bt="
    print "(d16.6)",bt

    

end program array1
