
module problem

    ! Specify the function f defining the problem to be solved, 
    ! and the true solution for checking the answers.

    ! These functions must be consistent: the second derivative of
    ! u_true should equal -f(x).


contains

    real(kind=8) function f(x)
        ! right hand side
        implicit none
        real(kind=8), intent(in) :: x
        f = 100.d0 * exp(x)
    end function f

    real(kind=8) function u_true(x)
        ! true solution
        implicit none
        real(kind=8), intent(in) :: x
        u_true = (100.d0*exp(1.d0) - 60.d0)*x + 120 - 100.d0*exp(x)
    end function u_true


end module problem
