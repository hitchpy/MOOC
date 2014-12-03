

.. _lab11:

Lab 11: Tuesday May 6, 2014
=============================

* Note that there are several example codes in the class repository that
  might be useful, e.g.

   * `$UWHPSC/codes/openmp`
   * `$UWHPSC/codes/lapack/random`
   * `$UWHPSC/2013/homeworks`
   * `$UWHPSC/2013/solutions`
    
* Discussion of random number generators.
  See `UWHPSC/labs/lab11`
* Questions about OpenMP?

Programming problem
---------------------------

Work on this in groups!

#. Write a Fortran program to do the following:

   * input `seed` and `n` from the command line
   * seed the random number generator using `init_random_seed` from the
     `random_util.90` module,
   * generate an array `x` of `n` random numbers
   * compute the mean of these values:  the sum of all elements of `x`
     divided by `n`.  Do this with a `do` loop.
   Since `random_number` produces numbers that should be uniformly
   distributed between 0 and 1, the mean should be approximately 0.5
   for large `n`.  It can also be shown that for a uniform
   distribution, the difference between the mean of a sample of `n` numbers
   and the true mean of the distribution should decay to zero like
   :math:`1/\sqrt{n}` as :math:`n\longrightarrow\infty`.  Do you observe this?

#. Modify your code to use OpenMP by using an `omp parallel do` loop
   with a suitable reduction to compute the mean.

#.  There is a quiz for this lab.
