
.. _homework4:

==========================================
Homework 4 [2014]
==========================================


Due Tuesday, May 27, 2014, by 11:00pm PDT.

**Part 1**

#. The directory `$UWHPSC/homeworks/homework4/part1` contains a possible
   solution to part of Homework 3.  The subroutine `approx_norm.f90`
   implements only OpenMP Method 1 from that homework.

   Modify `approx_norm.f90` to implement coarse grain parallelism in which
   the samples are manually split up into groups based on the number of threads
   and each thread handles  a set of `k` from its own
   `kstart` to `kend`, through the use of parallel blocks.  Each thread
   should  keep track of its on private `max_ratio_thread` for the maximum
   ratio it sees.
   Among all the threads, they should handle `k = 1` up to `nsamples`.
   After computing the maximum over its assigned set of `k` values,
   each thread should update a global `max_ratio` inside a critical block
   in order to avoid race conditions.

   You might want to base your code on `$UWHPSC/codes/openmp/pisum2.f90`,
   for example.

   You will need to know `nthreads`, the number of threads being used
   in the subroutine `approx_norm`.  Do not change the calling sequence of
   this subroutine or add additional module variables.  Instead call
   `omp_get_num_threads` at an appropriate point in the subroutine.

   You should not need to modify the main program, the `random_util`
   module, or the `Makefile`.  

   After modifying `approx_norm.f90`, it should then give output like this::

        $ make test
        gfortran -O2 -fopenmp -c  random_util.f90 
        gfortran -O2 -fopenmp -c  approx_norm.f90 
        gfortran -O2 -fopenmp -c  main.f90 
        gfortran -O2 -fopenmp random_util.o approx_norm.o main.o -o test.exe
        Wrote data to input_data.txt
        ./test.exe
         nthreads =            2
         seed1 for random number generator:       12345
        Thread  1 will take k =    501 through k =   1000
        Thread  0 will take k =      1 through k =    500
        CPU time =   0.00330300 seconds
        Elapsed time =   0.00209300 seconds
        Using matrix with n =    50   True 1-norm:          28.091553
        Based on   1000 samples, Approximate 1-norm:        24.998100

   Note that there are timing commands in the main program.  Depending
   on what system you are running on, you won't get the same times.
   You can experiment with whether using more threads improves the elapsed
   time, e.g. ::

        $ make test n=200 nsamples=10000 nthreads=1
        Wrote data to input_data.txt
        ./test.exe
         nthreads =            1
         seed1 for random number generator:       12345
        Thread  0 will take k =      1 through k =  10000
        CPU time =   0.42308200 seconds
        Elapsed time =   0.42331299 seconds
        Using matrix with n =   200   True 1-norm:         110.150308
        Based on  10000 samples, Approximate 1-norm:       100.499256


        $ make test n=200 nsamples=10000 nthreads=2
        Wrote data to input_data.txt
        ./test.exe
         nthreads =            2
         seed1 for random number generator:       12345
        Thread  1 will take k =   5001 through k =  10000
        Thread  0 will take k =      1 through k =   5000
        CPU time =   0.45311000 seconds
        Elapsed time =   0.25193700 seconds
        Using matrix with n =   200   True 1-norm:         110.150308
        Based on  10000 samples, Approximate 1-norm:       100.499256
    
  But see also the comments in Part 2 below about using OpenMP together
  with random number generators.

**Part 2**

#. Write a Fortran 90 code to solve the Gambler's Ruin problem considered
   in :ref:`lab13`, which is also visible at
   `here <http://nbviewer.ipython.org/url/faculty.washington.edu/rjl/notebooks/GamblersRuin.ipynb>`_.

   Base the code on the IPython notebook example function 
   `walk`, with one change: do not keep track of the path, so you do not
   need a parameter `n1path`.  

   A main program `main1.f90` can be found in
   `$UWHPSC/homeworks/homework4/part2`, along with `Makefile1`. 
   Write the module `gamblers.f90` that is needed, based on the template
   `$UWHPSC/homeworks/homework4/part2/gamblers.f90` -- do not change the
   module variables or the calling sequence of the code.  

   In the `walk` subroutine, you could either call `random_number(r)` for a
   scalar `r` on every step, or use an array of random numbers that is filled
   with `max_steps` values
   by one call to `random_number` before starting the loop, similar to the
   code provided for `approx_norm.f90` in Part 1.  
   **Use the second approach** with a single call to `random_number` because
   this makes the results more repeatable once OpenMP is added in the next
   part.  Threads might call
   `random_number` in a different order if you run the code multiple times,
   even with the same seed, and you want each walk to have the same set of
   random numbers.  They will if each calls `random_number` only once per
   walk, but not necessarily if each walk calls it multiple times.

   Which approach gives faster code depends on whether `max_steps`
   is reasonable for the problem being solved (in which case it can
   be faster to call it once than to have a subroutine call every
   step) or is much larger than needed (in which case it might be
   slower because many more random numbers are generated than needed).

   When you have filled in this module properly, you should get results like::

        $ make test1 -f Makefile1
        Wrote data to input_data.txt
        ./test1.exe
        n1 =      4
        n2 =      4
        p =   0.6000
        max_steps =      10000
         seed1 for random number generator:        1111
        In step      1 r =  0.7949 and n1 =      3 n2 =      5
        In step      2 r =  0.5090 and n1 =      4 n2 =      4
        In step      3 r =  0.2824 and n1 =      5 n2 =      3
        In step      4 r =  0.7906 and n1 =      4 n2 =      4
        In step      5 r =  0.7094 and n1 =      3 n2 =      5
        In step      6 r =  0.0509 and n1 =      4 n2 =      4
        In step      7 r =  0.1227 and n1 =      5 n2 =      3
        In step      8 r =  0.4534 and n1 =      6 n2 =      2
        In step      9 r =  0.2900 and n1 =      7 n2 =      1
        In step     10 r =  0.3408 and n1 =      8 n2 =      0
        Stopped after     10 steps with n1 =      8, n2 =      0
        After     10 steps, the winner is player  1

#. Add another main program `main2.f90` that uses an omp parallel do loop to 
   take many random walks and compute the fraction of wins by each player,
   and also the average number of steps in the walk.  You can use the Python
   code in `$UWHPSC/labs/lab13/GamblersRuin.ipynb` as a model for how to do
   this.  Keep the following in mind:
   
   * Within the loop you will call the `walk` function you wrote for
     the previous problem, and the `gamblers.f90` module should not have to
     change at all.  

   * Create a second Makefile named `Makefile2` with::

           OBJECTS = random_util.o gamblers.o main2.o

     and a target `test2` that runs the new version of the code.

     The Makefile should also set two additional parameters `kwalks` with default
     value 500 and `nthreads` with default value 2.  These values should be
     written to `input_data.txt` as part of the work done for the `data`
     target.  

     The main program should also print these values out.

   * Also keep track of the number of walk steps taken by each thread
     by introducing an array `nsteps_thread` and print these out at
     the end of the program.  

   * Add timing, using both `cpu_time` and `system_clock`, to time the
     main loop over `k = 1,kwalks`.  You can copy the necessary code
     from `part1/main.f90`.  (Don't forget to declare the necessary
     variables.)

   Sample output (with parameters giving longer walks)::

        $ make test2 -f Makefile2 n1=50 n2=50 p=0.501 nthreads=4
        Wrote data to input_data.txt
        ./test2.exe
        n1 =     50
        n2 =     50
        p =   0.5010
        kwalks =    500
        max_steps =     10000
        nthreads =  4
         seed1 for random number generator:        1111
        CPU time =   0.10481100 seconds
        Elapsed time =   0.08621634 seconds
        Warning:      3 walks out of      500 did not result in a win by either
        player
         
        Player 1 won  0.5600 fraction of the time, Player 2 won  0.4340 fraction of
        the time
        True probabilities are P1 =  0.5498 P2 =  0.4502
        The average path length is   2457
        True mean path length is   2491
         Total number of steps taken by each thread:
          Thread         0 took    312642 steps
          Thread         1 took    305352 steps
          Thread         2 took    290626 steps
          Thread         3 took    320061 steps

   **Note on steps per thread:** Even if you make sure each walk gets
   the same set of random numbers by calling `random_number` only once
   per walk, the threads might split up the walks differently if you
   run the code repeatedly, so the last set of numbers above could change
   but the computed fraction of wins and average path length should not
   change if you keep running with the same seed.

   **Note on timings:**  With this code you will probably not see any
   speed up due to the use of OpenMP even if it appears the work is
   evenly divided.  This is because each call to `walk` by either thread
   requires a call to `random_number` and the random number generator is
   thread safe, meaning that it contains critical blocks so that only one
   thread at a time can be accessing its internal state.  This makes the
   code usable with OpenMP, but if much of the work being done is in
   the random number generator (as in this simple code) then it may not run
   much faster and perhaps even slower than just using one thread.

   A possible way to speed up the code would be to generate all the random
   numbers needed for all the walks in the main program before the loop on
   `kwalks`, but this would require more re-writing of code.  You are
   welcome to experiment with this if you wish, but not required.
   Turn in code that follows the instructions above for the assignment.


To submit
---------

* At the end, you should have committed the following 
  files to your repository:

  **Part 1**

  * `$MYHPSC/homework4/part1/Makefile`  (unchanged from original)
  * `$MYHPSC/homework4/part1/main.f90`  (unchanged from original)
  * `$MYHPSC/homework4/part1/random_util.f90`  (unchanged from original)
  * `$MYHPSC/homework4/part1/approx_norm.f90`

  **Part 2**

  * `$MYHPSC/homework4/part2/Makefile1`  (unchanged from original)
  * `$MYHPSC/homework4/part2/random_util.f90`  (unchanged from original)
  * `$MYHPSC/homework4/part2/main1.f90`  (unchanged from original)
  * `$MYHPSC/homework4/part2/gamblers.f90`
  * `$MYHPSC/homework4/part2/main2.f90`
  * `$MYHPSC/homework4/part2/Makefile2`

  Note that we should be able to run your code by giving commands like
  those given above.  But also if we write a new main program that calls 
  your subroutine `approx_norm` or `walk`, that should also work.

  **Please be sure you have the specified directory and file names.**
  It is hard to grade otherwise, and points will be deducted.
  

  Make sure you push to bitbucket after committing.

* Submit the commit number that you want graded by following the link
  provided on the `Canvas page for Homework 4
  <https://canvas.uw.edu/courses/893991/assignments/2513531>`_.


