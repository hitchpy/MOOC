subroutine area(r,a)

    use constants, only: pi
    implicit none
    real(kind=8), intent(in) :: r
    real(kind=8), intent(in) :: a


    a = pi * r**2
    
end subroutine area
