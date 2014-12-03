
.. _homework3:

==========================================
Homework 3 [2014]
==========================================


Due Thursday, May 8, 2014, by 11:00pm PDT.

Recall from linear algebra that if :math:`\|x\|` is some vector norm then
the corresponding matrix norm is defined by

    :math:`\|A\| = \max_{x\neq 0} \frac{\|Ax\|}{\|x\|}`

This means that :math:`\|Ax\| \leq \|A\|\|x\|` for all :math:`x`, and in fact
the norm of :math:`A` is the smallest value :math:`c\geq 0` such that the
inequality :math:`\|Ax\| \leq c\|x\|` holds for *all* vectors :math:`x` of the
appropriate shape.

Recall also that the 1-norm of a vector of length :math:`n` is defined by

    :math:`\|x\|_1 = \sum_{i=1}^n |x_i|`

It is fairly easy to show that there is a simple formula for the
corresponding matrix 1-norm of an :math:`n \times n` matrix :math:`A`:

    :math:`\|A\|_1 = \max_{j=1,2,\ldots,n} \sum_{i=1}^n |A_{ij}|`

In other words, compute the 1-norm of each column vector and the 1-norm of
the matrix is the maximum of these column sums.

The goal of this homework is to:

* Compute the 1-norm by this formula and
  compare it to an approximation obtained by computing the maximum of
  :math:`\frac{\|Ax\|}{\|x\|}` over a large number of random vectors.

* Speed the code up using OpenMP, exploring two different ways to add
  parallelization to the approximation.


The directory `$UWHPSC/homeworks/homework3` contains some code to get you
started on this assignment:

* `main.f90` is a main program that reads in some data values,
  creates a random matrix, computes the 1-norm by the formula, and then
  calls `approx_norm1` and prints out the approximation.

* `approx_norm.f90` has a template for how this module should look, but does 
  not contain the necessary code.

* `random_util.f90` that has a function to seed the Fortran random number
  generator based on a single integer that you can set to some positive
  integer (so the results are repeatable if you run the code twice).
  If you set `seed = 0` then it will generate a "random" seed.
  See :ref:`random` for some references.

* `Makefile` is set up so to create a file `input_data.txt` that is read in
  by the Fortran code, to illustrate one way of changing parameters easily
  at the command line.  

For example::

    $ make test

will compile the code and run it with some default parameter values set in
the Makefile.  But you can also do, e.g.::

    $ make test n=200

to try a :math:`200 \times 200` random matrix instead of the default
:math:`50 \times 50`, or::

    $ make test n=200 seed=0

to also change the seed for the random number generator.  
(Also note that `seed=0` indicates
that it should give different random results every time you run this.)
    
Assignment
----------

#. Write the Fortran code needed in `approx_norm.f90` so that 
   the program gives sensible output, e.g.::

        $ make test

        Wrote data to input_data.txt
        ./test.exe
         seed1 for random number generator:       12345
        Using matrix with n =    50   True 1-norm:          28.091553
        Based on   1000 samples, Approximate 1-norm:        24.998100

   Note that the approximation to the 1-norm should always be less than
   or equal to the true 1-norm, since in the program you are only computing
   the maximum over some finite number of random vectors rather than over all
   :math:`x`.  (And in general this does **not** give a very good
   approximation, even when the number of samples is very large!)

   To write this code, note that you want to loop over `nsamples` different
   random vectors :math:`x`, and for each vector compute
   :math:`\|Ax\|/\|x\|`.  Keep track of the maximum of these as you go along.
   Note that `nsamples` is also read in from `input_data.txt` and a 
   default value is specified in the Makefile.

   You will need to determine the size of the matrix `a` passed into to the
   subroutine since that is not an explicit argument.  For this you can
   use::

        n = size(a,1)   ! number of rows in a


#. Change the `Makefile` so that it compiles with the `-fopenmp` flag and
   add OpenMP directives to the main program so that the loop on `j` for
   computing the true 1-norm is a parallel do loop.  Think about what variables
   need to be declared `private`.  Make sure the program still runs and
   gives consistent results whether you compile with OpenMP or not.

   **Note:**  To recompile with or without OpenMP, you should first do::

        $ make clean

   in the directory to make sure it recompiles all the `*.o` object files
   with the correct compiler flag.

#. Add an input parameter `nthreads` that is read in from `input_data.txt`
   and used to set the number of threads to use in the main program.  Also
   add this to the Makefile with a default value 2 and print it out from the
   main program, so that for example::

       $ make test nthreads=4
       Wrote data to input_data.txt
       ./test.exe
        nthreads =            4
        seed1 for random number generator:       12345
       Using matrix with n =    50   True 1-norm:          28.091553
       Based on   1000 samples, Approximate 1-norm:        24.998100

    
#. Parallelize `approx_norm1` using OpenMP.  
   Note that there is more than one way to do this.

   1. You could parallelize the outermost loop over random vectors 
   (so different threads are working on different vectors `x`), or 

   2. You could loop over the different :math:`x` vectors as in the
   serial code, but then parallelize the work that must be done in computing 
   :math:`\|Ax\|/\|x\|` for each `x`.  

   Implement both these approaches, and add a parameter `method` so that
   `method=1` means the first approach and `method=2` means the second
   approach.  Handle this parameter similar to the other input data, with
   a default value in the Makefile, and with `main.f90` reading it in from
   the file `input_data.txt`.

   Add `method` as a module variable to `approx_norm.f90` in order to pass
   the value from the main routine into the subroutine.  Do not change the
   calling sequence of the subroutine.

#. (You don't need to turn anything in for this part since timing parallel
   codes can be dicey on some machines.)

   Experiment with the two methods implemented above to see which approach
   seems to be better on large problems.   For example you might try::

        $ time make test n=50 nsamples=100000 nthreads=1 method=1
   
   and then see what happens as you increase the number of threads with this
   method, and then repeat with `method=2`.

   Note that this problem has small matrices and vectors but lots of samples.

   Also see what happens if the matrix is big but the number of samples is 
   relatively small, e.g. ::

        $ time make test n=5000 nsamples=100 nthreads=1 method=1

   Can you understand the behavior you see?  
   If you get counter-intuitive results, try to understand why.

To submit
---------

* At the end, you should have committed the following 
  files to your repository:

  * `$MYHPSC/homework3/Makefile`
  * `$MYHPSC/homework3/main.f90`
  * `$MYHPSC/homework3/random_util.f90`  (unchanged from original)
  * `$MYHPSC/homework3/approx_norm.f90`

  Note that we should be able to run your code by giving commands like
  those given above.  But also if we write a new main program that calls 
  your subroutine `approx_norm1`, that should also work.

  Make sure you push to bitbucket after committing.

* Submit the commit number that you want graded by following the link
  provided on the `Canvas page for Homework 3
  <https://canvas.uw.edu/courses/893991/assignments/2504886>`_.

