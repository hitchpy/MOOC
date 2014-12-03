
module constants

    ! define various constants
    implicit none
    real(kind=8) :: pi, e
    save

contains

    subroutine set_constants()
        implicit none
        pi = acos(-1.d0)
        e = exp(1.d0)
    end subroutine set_constants

end module constants
