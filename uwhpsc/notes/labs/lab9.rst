

.. _lab9:

Lab 9: Tuesday April 29, 2014
=============================

Programming problem
---------------------------

Work on this in groups!

Write a Fortran program to compute the roots of a quadratic equation and
determine the absolute and relative error in the roots computed.  The
program should do the following:

* Prompt the user and then read in the desired roots `x1true` and `x2true`
  (See the example below).
* Determine the coefficients `a, b, c` so that the quadratic equation
  :math:`a x^2 + bx + c =0` has the desired roots (you can set `a=1`).
* Using `a,b,c`, compute the roots `x1` and `x2` by using the quadratic 
  formula.
* Print out the "true" and computed values for each root along with the
  absolute and relative error in each.

So you should be able to do something like this::

    $ gfortran quadratic.f90
    $ ./a.out
     input x1true, x2true: 
    2.5, 3

    Coefficients:  a =     0.100000E+01  b =    -0.550000E+01  c = 0.750000E+01

    Root x1  computed:  0.250000000000000E+01  true:  0.250000000000000E+01
             absolute error:     0.000000E+00  relative error:     0.000000E+00
      
    Root x2  computed:  0.300000000000000E+01  true:  0.300000000000000E+01
             absolute error:     0.000000E+00  relative error:     0.000000E+00

Don't worry too much about the formatting but you might want to print out 15
digits in the computed roots.

You might want to assume the values are entered with `x1true <= x2true` so
you know which root from the quadratic equations goes with which original
value. (And print out an informative error message otherwise.)

Test it out
-----------

Test a variety of values to see that it's working.  

Once it's working on reasonable values, try the following: 

* `x1 = 1e-12,  x2 = 2`  
* `x1 = -2, x2 = 1e-12`

In each case you should find that one root is computed accurately
but the other root has a large relative error (few digits of accuracy).

Figure out why "catastrophic cancellation" is the problem.

Improvements
------------

* Improve the code by noticing that if one root is calculated accurately, the
  other root can be calculated from the fact that `x1 * x2 = c`.

* Remove the assumption that `x1true <= x2true`.


