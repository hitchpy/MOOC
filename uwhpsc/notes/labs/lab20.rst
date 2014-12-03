

.. _lab20:

Lab 20: Thursday June 5, 2014
=============================


* The directory `$UWHPSC/labs/lab20` contains a Fortran code
  `fourier_sum.f90` that computes terms in the Fourier sine series
  for the function :math:`f(x) = e^x` and prints them out to the file
  `frames.txt` along with partial sums of the series.  
  It also creates a second file `xf.txt` that contains the values of `x`
  and `f(x)` on the finite grid where the terms and sum are computed.

  The Python script `animate.py` reads in these two files and produces
  an animation, which can be viewed at
  `<http://faculty.washington.edu/rjl/FourierSum.html>`_.

  We will go through these codes.

* You may have to install JSAnimation on your SMC project, see
  :ref:`animation`.


   
Work on in pairs:
------------------

* The subdirectory `$UWHPSC/labs/lab20/gamblers_ruin` contains a modified
  version of the solution to Homework 4, part 2,  and does a single random
  walk.  The `gamblers.f90` code has been modified to print `n1, n2` each
  step to the file `walk.txt`.

  Produce a script `animate.py` and modify `Makefile1` so that you can do::
  
        $ make movie -f Makefile1 n1=10 n2=10 seed=1234

  and generate a movie like the one shown at
  `<http://faculty.washington.edu/rjl/RandomWalk.html>`_.

* Don't worry about the Makefile at first, get `animate.py` working.

* For testing, you might want to use smaller `n1` and `n2` so there are 
  fewer `png` files to generate.


**There is quiz for Lab 20**


