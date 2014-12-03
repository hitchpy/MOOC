program demo

    use constants, only: set_constants, pi, e

    implicit none
    real(kind=8) :: r, a

    call set_constants()
    print *, "In main program, pi = ", pi, "  e = ", e
    print *, "Input r: "
    read *, r
    call area(r,a)
 10 format("A circle of radius ", e21.15," has area ", e21.15)
    print 10, r, a

end program demo
