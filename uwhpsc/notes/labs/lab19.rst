

.. _lab19:

Lab 19: Tuesday June 3, 2014
=============================


* We will go through the notebook `$UWHPSC/homeworks/project/Heat_Equation.ipynb`, also 
  visible at
  `<http://nbviewer.ipython.org/url/faculty.washington.edu/rjl/notebooks/Heat_Equation.ipynb>`_.

  This notebook gives a brief introduction to the heat equation and
  two numerical methods for its solution, an explicit method and the
  more stable implicit Crank-Nicolson method.
   
Some things to try
------------------

* You might want to make a copy of this notebook before you start playing
  with it.

* Experiment with different initial conditions for the heat equation.

* Create an animation (in the notebook) of the numerical solution to the 
  heat equation along with the true solution.

* The Crank-Nicolson method is *second order accurate*: the error should
  go to zero like :math:`{\cal O}(\Delta t^2 + \Delta x^2)` as the grid is
  refined.  So increasing both `n` (the number of spatial points) and
  `nsteps` (the number of time steps) by a factor of 2 should reduce the
  error by a factor of 4.  Test this out.

* Compute or look up the Fourier sine series for some interesting function
  and try this as initial conditions for the heat equation.  Compare the
  true solution with the numerical solution (where the "true solution" might
  be estimated by adding up a large but finite number of terms in the 
  Fourier series).

* Try using Sympy to compute the coefficients in the Fourier sine series.

**There is quiz for Lab 19**


