
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
    real(kind=8), allocatable :: r(:)
    integer :: nstep, n1, n2

    allocate(r(max_steps))
    call random_number(r)

    ! initialize n1 and n2 with input values
    ! need to change n1 and n2 internally during this walk but do not
    ! want to affect n1 and n2 in main program.
    n1 = n1in
    n2 = n2in

    open(unit=22, file='walk.txt', status='unknown')
    write(22,222) n1,n2
222 format(2i6)

    do nstep=1,max_steps
        if (r(nstep) <= p) then
            n1 = n1 + 1
            n2 = n2 - 1
        else
            n1 = n1 - 1
            n2 = n2 + 1
        endif 
        write(22,222) n1,n2
        if (verbose) then
            print 601, nstep,r(nstep),n1,n2
601         format("In step ",i6, " r = ",f7.4, " and n1 = ",i6, " n2 = ",i6)
            endif
        if (min(n1,n2) == 0) then
            exit
            endif

        enddo

    nsteps = nstep  ! to return number of steps taken

    if (verbose) then
        print 602, nstep,n1,n2
602     format("Stopped after ",i6," steps with n1 = ",i6,", n2 = ",i6)
        endif

    if (n1 == 0) then
        winner = 2
    else if (n2 == 0) then
        winner = 1
    else
        winner = 0
    endif 

    end subroutine walk

end module gamblers
