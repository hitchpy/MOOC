

.. _lab16:

Lab 16: Thursday May 22, 2014
=============================


Problem to solve
----------------

    
* Adapt the program `$UWHPSC/codes/lapack/randomsys3.f90` to use a specific matrix A in place of
  the random matrix used in the original code.  The matrix to use is the Hilbert matrix defined by

    :math:`a_{i,j} = \frac{1}{i+j-1}`

  This is a notorious matrix since it is always nonsingular but is very ill-conditioned even for
  moderately small values of `n`.

  For more discussion of this matrix, and a formula for how the condition number grows with `n`,
  see this `Cleve's Corner blog post
  <http://blogs.mathworks.com/cleve/2013/02/02/hilbert-matrices/#73083b00-1b97-4570-a516-31796a031dc4>`_.

* Note that in order to create an executable for your program, in the linking step you will need
  to make sure `gfortran` also links in the BLAS and LAPACK library.  See the `LFLAGS` set in
  `$UWHPSC/codes/lapack/Makefile for the arguments you need to add to the linking step.

* Instead of using the `random_number` subroutine to generate a random `x` for checking the
  relative error, as is done in `$UWHPSC/codes/lapack/randomsys3.f90`, try taking `x` to be 
  a vector of all 1's.  (And as in the original code compute :math:`b = Ax` usint `matmul` and
  then solve the system to recover `x`.)  Print out the computed `x` as well as computing the
  relative error in the 1-norm as in the original code.  How well does it do?  How does the
  accuracy relate to the condition number?

* You might want to look at the `dgecon documentation
  <http://www.netlib.no/netlib/lapack/double/dgecon.f>`_.

* Try different values of `n` with your program to see if it gives the expected behavior.
  Note that the LAPACK function `dgecon` does not compute the exact condition number but only
  estimates it.  Also note that the program estimates the 1-norm condition  number, while the
  approximate formula is for the 2-norm condition number (but they grow in a similar exponential 
  fashion).

**If you have time to do more...**

* Modify your code by creating a Fortran function `hilbert_condition` that returns the condition number
  estimate for a given value of `n`.  

  Then write a main program that loops over `n` from 1 to 20, computes the estimate for each `n`, 
  and writes a text file with two columns `n` and the estimate.  These statements might be
  useful::

            open(21, file='cond.txt',status='unknown')

            do n=1,20
                cond = hilbert_condition(n)
                ! print *, "cond = ",cond
                write(21, 210) n,cond
         210    format(i4,e16.6)
                enddo


* The text file produced should be readable by the Python script 
  `$UWHPSC/labs/lab16/plot_cond.py`, which plots the results on a logarithmic scale, along with
  what the formula predicts.

* For the function version you do not need to solve a linear system, so you don't need to call
  `dgesv`, but you do need to compute the LU factorization of A before calling `dgecon`.   
  The could be done by calling `dgetrf` instead of `dgesv`.
  You might want to look at the `dgetrf documentation
  <http://www.netlib.no/netlib/lapack/double/dgetrf.f>`_.


**There is quiz for Lab 16**
