

.. _lab10:

Lab 10: Tuesday May 1, 2014
=============================

Programming problem
---------------------------

Work on this in groups!

#. The OpenMP code `$UWHPSC/labs/lab10/array_omp.f90` contains some bugs.
   Find the bugs and fix them so that it runs and gives output like this::

        $ gfortran -fopenmp array_omp.f90
        $ ./a.out
         nthreads =            6
         b and bt should be equal 
         b=
            0.270000D+02
            0.330000D+02
            0.390000D+02
            0.450000D+02
            0.510000D+02
         bt=
            0.270000D+02
            0.330000D+02
            0.390000D+02
            0.450000D+02
            0.510000D+02

#.  If :math:`A` is an :math:`n \times n` matrix and :math:`x` is a vector of
    length :math:`n`, then :math:`x^TAx` is a scalar, a "quadratic form"
    since it is the sum of terms of the form :math:`a_{ij}x_ix_j` that are
    quadratic in the elements of :math:`x` .

    Write an OpenMP code to compute this for a given matrix and vector.  Write
    out the matrix-vector multiplies as loops and use "omp parallel do" loops to
    compute first the vector :math:`Ax` and then the inner product of this with
    the vector `x`.  Test your code using the :math:`10 \times 10` identity matrix 
    for :math:`A` and :math:`x_i = i`, in which case the correct answer can be
    determined to be 385 from the formula

        :math:`\sum_{i=1}^n i^2 = \frac{n(n+1)(2n+1)}6.`


#.  There is a quiz for this lab.
