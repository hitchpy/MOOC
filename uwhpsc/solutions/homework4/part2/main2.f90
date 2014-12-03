
program main2

    use omp_lib
    use random_util, only: init_random_seed
    use gamblers, only: kwalks, max_steps, walk

    implicit none
    real(kind=8) :: p, q, p1, p2, frac1wins, frac2wins
    integer :: wins(0:2)
    integer :: seed, i, j, nthreads, k, n1, n2, nsteps, winner, thread_num
    integer :: nsteps_total 
    real(kind=8) :: nsteps_average, mean_length
    integer, allocatable :: nsteps_thread(:)
    integer(kind=8) :: tclock1, tclock2, clock_rate
    real(kind=8) :: t1, t2, elapsed_time


    open(unit=21, file='input_data.txt', status='old')
    read(21,*) n1
    read(21,*) n2
    read(21,*) p
    read(21,*) max_steps
    read(21,*) seed
    read(21,*) kwalks
    read(21,*) nthreads

    print "('n1 = ',i6)", n1
    print "('n2 = ',i6)", n2
    print "('p = ',f8.4)", p
    print "('kwalks = ',i6)", kwalks
    print "('max_steps = ',i9)", max_steps

    allocate(nsteps_thread(0:nthreads-1))
    nsteps_thread = 0

    !$ call omp_set_num_threads(nthreads)
    print "('nthreads = ',i2)", nthreads

    call init_random_seed(seed)

    wins = 0
    nsteps_total = 0

    !---- check the time for the main loop:

    call system_clock(tclock1)
    call cpu_time(t1)   ! start cpu timer

    !$omp parallel do private(k,nsteps,winner,thread_num) &
!   !$omp             reduction(+:nsteps_total) schedule(dynamic,1)
    !$omp             reduction(+:nsteps_total) 
    do k=1,kwalks
        call walk(n1, n2, p, .false., nsteps, winner)
        wins(winner) = wins(winner) + 1
        !$ thread_num = omp_get_thread_num()
        nsteps_thread(thread_num) = nsteps_thread(thread_num) + nsteps
        nsteps_total = nsteps_total + nsteps
        enddo

    call cpu_time(t2)   ! end cpu timer
    print 10, t2-t1
 10 format("CPU time = ",f12.8, " seconds")

    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")

    if (wins(0) > 0) then
        print 601, wins(0), kwalks
 601    format("Warning: ",i6, " walks out of ",i8, &
               " did not result in a win by either player")
        endif

    frac1wins = wins(1) / float(kwalks)
    frac2wins = wins(2) / float(kwalks)
    print 602, frac1wins,frac2wins
602 format(/,"Player 1 won ", f7.4," fraction of the time, Player 2 won ", &
           f7.4, " fraction of the time")
            
    q = 1.d0 - p
    p1 = (1.d0-(q/p)**n1) / (1.d0-(q/p)**(n1+n2))
    p2 = 1.d0 - p1
    print 603, p1, p2
603 format("True probabilities are P1 = ", f7.4," P2 = ",f7.4)

    nsteps_average = nsteps_total / float(kwalks)
    print 604, nsteps_average
604 format("The average path length is ",f9.2)

    mean_length = (n1/(q-p)) - ((n1+n2)/(q-p))*(1-(q/p)**n1)/(1-(q/p)**(n1+n2))
    print 605, mean_length
605 format("True mean path length is ",f9.2)

    print *, "Total number of steps taken by each thread:"
    do thread_num = 0,nthreads-1
        print 606, thread_num, nsteps_thread(thread_num)
606     format("  Thread ",i9," took ",i9," steps")
        enddo

end program main2
