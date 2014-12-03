
module bvp_solvers

contains

subroutine solve_bvp_direct(x, u_left, u_right, u)
    use problem, only: f
    implicit none
    real(kind=8), intent(in) :: x(0:)
    real(kind=8), intent(in) :: u_left, u_right
    real(kind=8), intent(out) :: u(0:)

    ! decleare local variables and continue writing this code...

end subroutine solve_bvp_direct

end module bvp_solvers
