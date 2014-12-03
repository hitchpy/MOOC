

.. _lab12:

Lab 12: Thursday May 8, 2014
=============================


Programming problem
---------------------------

Work on this in groups!

#. In Lab 11 you worked on a program to compute the mean of n random
   numbers.  A sample solution can be found at `$UWHPSC/labs/lab11/mean.f90`.

   Write a Fortran program that runs over different values of `n`, 
   and for each `n` generates a vector `x` containing `n` random numbers
   and then computes the mean of these.  Also compute the fraction of the
   numbers that lie in the first quartile (the fraction of `x(i)` values
   that are between 0 and 0.25) and the fraction that lie in the fourth
   quartile (between 0.75 and 1.0).  Since the `random_number` routine
   returns numbers uniformly distributed between 0 and 1, we expect  each of
   these fractions to be about 0.25.

   Use OpenMP to make the loop on `i` from 1 to `n` into a parallel do loop.

   Running this code should give something like this if you take as the `n`
   values :math:`n = 10^k` for :math:`k=2,3,\ldots,8`::


         Number of threads:            2
         input seed
        12345
         seed1 for random number generator:       12345

               n          mean        quartile 1    quartile 4
               100     0.51902466     0.22000000     0.24000000
              1000     0.47476778     0.27800000     0.22500000
             10000     0.49606601     0.25670000     0.25190000
            100000     0.50121669     0.24815000     0.25130000
           1000000     0.50001034     0.24986300     0.24979800
          10000000     0.49998532     0.24994350     0.24992770
         100000000     0.49995944     0.25003764     0.24995608


#. If you haven't already, study the code in
   `$UWHPSC/codes/openmp/pisum2.f90`
   and make sure you understand how this coarse grain parallelism works.
   Discuss with your neighbors.

#. If you have time, try to follow this model to make your code that
   computes the mean and quartiles work in a similar manner, where you
   break up the different values of `n` to be tested between different
   threads, e.g. in the above example one thread would take the
   first three values of `n` and the second thread would take the final
   two values of `n`.

#. Discuss with your neighbors whether this is a sensible way to try
   to use two threads on this problem.

#. There is a quiz for this lab.

#. Sample solutions can now be found in `$UWHPSC/labs/lab12`.
