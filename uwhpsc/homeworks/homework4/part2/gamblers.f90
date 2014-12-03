
module gamblers

    implicit none
    integer :: kwalks, max_steps
    save

contains

    subroutine walk(n1in, n2in, p, verbose, nsteps, winner)

    implicit none
    integer, intent(in) :: n1in,n2in
    real(kind=8), intent(in) :: p
    logical, intent(in) :: verbose
    integer, intent(out) :: nsteps, winner

    ! local variables
    real(kind=8) :: r
    integer :: nstep, n1, n2

    ! initialize n1 and n2 with input values
    ! need to change n1 and n2 internally during this walk but do not
    ! want to affect n1 and n2 in main program. 
    n1 = n1in
    n2 = n2in

    ! CONTINUE WRITING THIS SUBROUTINE....

    end subroutine walk

end module gamblers
