

.. _lab18:

Lab 18: Thursday May 29, 2014
=============================


* We will go through the notebook `$UWHPSC/homeworks/project/BVP.ipynb`, also 
  visible at
  `<http://nbviewer.ipython.org/url/faculty.washington.edu/rjl/notebooks/BVP.ipynb>`_.
  This outlines a recursive `domain decomposition
  <https://www.google.com/search?q=domain+decomposition&rlz=1C5CHFA_enUS534US534&espv=2&source=lnms&tbm=isch&sa=X&ei=R4GHU4uFKI-XyAT4t4Bo&ved=0CAYQ_AUoAQ&biw=1440&bih=779>`_  approach to solving a
  boundary value problem.  Part 1 of the project is to convert this into
  Fortran with OpenMP.
   
* Working in pairs, copy this notebook to `BVP2.ipynb` and modify it to solve
  a `Helmholtz equation <http://en.wikipedia.org/wiki/Helmholtz_equation>`_
  (in one dimension) of the form:

    :math:`u''(x) + k^2 u(x) = -f(x)`

  on the interval :math:`0<x<1` with specified boundary conditions.  

  As an exact solution, consider the case :math:`f(x)=0` in which case
  the general solution to :math:`u''(x) = -k^2 u(x)` is 
  :math:`u(x) = c_1 \cos(kx) + c_2 \sin(kx)`.

  The boundary value problem has a unique exact solution for any boundary
  values :math:`u(0)` and :math:`u(1)` provided that :math:`k` is not an 
  integer multiple of :math:`\pi`.  (Insert :math:`x=0` and :math:`x=1` into
  the general solution and determine :math:`c_1` and :math:`c_2` so that the
  boundary conditions are satisfied.)

  You might try values such as::

        k = 15.
        u_left = 2.
        u_right = -1.

  You will need to use at least 40 or 50 grid points to get a solution that
  looks at all reasonable.
  If you make :math:`k` larger, the solution will be more oscillatory and
  you will need even more grid points to get a reasonable approximation.

* Work through as much of the notebook as you can, adjusting things to
  solve the Helmholtz equation.  The main objective is to work through the
  notebook and understand what is being done.

  Some tips:

  * Add another parameter `k` to the `solve_BVP_*` functions.  

  * In setting up the tridiagonal matrix in `solve_BVP_direct`, you will need
    to modify the diagonal terms for the difference equation that
    approximates the Helmholtz equation,

    :math:`\frac{U_{i-1} - 2U_i + U_{i+1}}{\Delta x^2} + k^2 U_i = -f(x_i)`

    This gives the linear system to be solved:

    :math:`U_{i-1} + (k^2\Delta x^2 - 2) U_i + U_{i+1} = -\Delta x^2 f(x_i)`

    along with the boundary conditions.

  * If you get to the divide-and-conquer approach, you will have to modify
    the equation for the mismatch to take into account the modification to
    the linear system being solved.

* There is now a sample solution in the repository, visible at
  `<http://nbviewer.ipython.org/url/faculty.washington.edu/rjl/notebooks/BVP_helmholtz.ipynb>`_.

* **There is quiz for Lab 18**


