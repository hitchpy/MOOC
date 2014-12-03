
module problem

    ! Save the value of k and define the true solution.

    integer :: k
    save

contains

    real(kind=8) function u_true(x,t)
        ! true solution
        implicit none
        real(kind=8), intent(in) :: x, t
        u_true = exp(-k**2 * t) * sin(k*x)
    end function u_true


end module problem
