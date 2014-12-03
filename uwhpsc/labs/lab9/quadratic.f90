
! This version illustrates catastrophic cancellation if you give it input
! such as x1true = 1e-12, x2true = 2.

program quadratic

    implicit none
    real(kind=8) :: a,b,c,x1,x2,x1true, x2true, err1, err2, s
    print *, "input x1true, x2true: "
    read *, x1true, x2true

    if (x1true > x2true) then
        print *, "Require x1true <= x2true"
        stop
        endif

    a = 1.d0
    b = -(x1true + x2true)
    c = x1true*x2true

	print 600, a,b,c
600 format(/,"Coefficients:  a = ",e16.6, "  b = ",e16.6, "  c = ",e16.6,/)

    s = sqrt(b**2 - 4.d0*a*c)
    x1 = (-b - s) / (2.d0*a)
    x2 = (-b + s) / (2.d0*a)

    err1 = x1 - x1true
    err2 = x2 - x2true

    print 601, 1, x1, x1true
    print 602, abs(err1), abs(err1/x1true)
    print *, ' '
    print 601, 2, x2, x2true
    print 602, abs(err2), abs(err2/x2true)

601 format("Root x",i1,"  computed: ",e22.15,"  true: ",e22.15)
602 format("         absolute error: ",e16.6,"  relative error: ",e16.6)


end program quadratic
