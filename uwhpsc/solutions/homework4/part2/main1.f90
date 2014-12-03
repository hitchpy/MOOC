
program main1

    use random_util, only: init_random_seed
    use gamblers, only: max_steps, walk

    implicit none
    real(kind=8) :: p, q, p1, p2, frac1wins, frac2wins
    integer :: seed, i, j, nthreads, k, n1, n2, nsteps, winner


    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n1
    read(21,*) n2
    read(21,*) p
    read(21,*) max_steps
    read(21,*) seed

    print "('n1 = ',i6)", n1
    print "('n2 = ',i6)", n2
    print "('p = ',f8.4)", p
    print "('max_steps = ',i9)", max_steps


    call init_random_seed(seed)

    call walk(n1, n2, p, .true., nsteps, winner)

    if (winner==0) then
        print 601, max_steps
601     format("There is no winner after ",i6," steps")
    else
        print 602, nsteps, winner
602     format("After ",i6," steps, the winner is player ",i2)
        endif

end program main1
