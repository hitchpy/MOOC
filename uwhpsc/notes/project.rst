
.. _project:

==========================================
Final Project [2014]
==========================================

Due Wednesday, June 11, 2014, by 11:00pm PDT.

.. warning:: A survey link has been added at the end of this page.

**Some resources and hints:**  See :ref:`project_hints`.  

**Part 1**

The notebook `$UWHPSC/project/BVP.ipynb` illustrates
the solution of the two-point boundary value problem :math:`u''(x)
= -f(x)` by setting up and solving a tridiagonal system.  
This was discussed in :ref:`lab18`.
(`nbviewer
<http://nbviewer.ipython.org/url/faculty.washington.edu/rjl/notebooks/BVP.ipynb>`_, `video <https://canvas.uw.edu/courses/893991/wiki/lab-18>`_)

The goal of Part 1 is to create a Fortran 90 module that contains 
subroutines `solve_bvp_direct` and `solve_bvp_split` that mimic
the Python codes, and then add OpenMP to the latter in order to 
solve the four boundary value problems simultaneously with 4 threads.

You do **not** need to convert the recursive code to Fortran, although
if you want more of a challenge you might try to do this with OpenMP.

The directory `$UWHPSC/project/part1` contains the following files:

* `main1.f90` a driver program for the first part below.
* `problem.f90`, a module containing the problem description (right hand
  side function `f` and true solution function `u_true`).
* `bvp_solvers.f90`, the template for a module that should be filled in 
  properly with three subroutines following the directions below.
* `Makefile1` for the first part below.

Note that the problem specified in `problem.f90` is the same one used in the
IPython notebook and in Lecture 23.  You might want to test your routines
with a different problem, e.g. choose a different true solution and
differentiate twice to get the corresponding right hand side (and remember
to change the boundary conditions specified in the main program to be
suitable).  We might do the same to test your routines.

#.  Finish writing the subroutine `solve_BVP_direct`  in the module
    `bvp_solvers.f90`.  This routine should set up a
    tridiagonal system and solve it using the LAPACK subroutine `DGTSV`.
    As input it takes a 1-dimensional array `x`, and two real numbers
    `u_left` and `u_right` and the output should be an array `u` containing
    the approximate solution at the points specified in `x`.  It can be assumed
    that the points in `x` are equally spaced.

    **Notes:**

    * For consistency with the Python code, the arrays `x` and `u`
      have been declared to start with index 0 rather than 1.  Note how this
      is indicated in the declaration of these arrays in the subroutine
      `solve_bvp_direct`.
    
    * Instead of `DGTSV <http://www.netlib.no/netlib/lapack/double/dgtsv.f>`_,
      you could use 
      `DPTSV <http://www.netlib.no/netlib/lapack/double/dptsv.f>`_, 
      which solves a tridiagonal
      `symmetric positive definite
      <http://en.wikipedia.org/wiki/Positive-definite_matrix>`_ 
      system.  The matrix `A` for this boundary
      value problem is symmetric but is negative definite rather than
      positive definite, so if you want to do this you would have to negate
      the system and solve :math:`-u''(x) = f(x)`.  
      The `DPTSV` routine is slightly more efficient, but it's fine to
      use `DGTSV`.

    **NOTE: The output below was created using `DPTSV` and if you use 
    `DGTSV` slightly different results are obtained.**

    Add a print statement to your subroutine so that running the code
    produces output like this::

        $ make test -f Makefile1 
        Wrote data to input_data.txt
        ./test.exe
         n =           20
        Solving tridiagonal system with n =     20
         error_max =   3.99812551471256938E-003


        $ make test -f Makefile1 n=10000
        Wrote data to input_data.txt
        ./test.exe
         n =        10000
        Solving tridiagonal system with n =  10000
         error_max =   1.69490235180091986E-008

#.  Add another subroutine `solve_bvp_split` to the `bvp_solvers` module
    that has the same calling sequence and produces the same result, but  
    that splits the interval roughly in half and makes four calls to
    `solve_bvp_direct`, following the same idea as in the IPython notebook.

    Provide a second main program `main2.f90` and makefile `Makefile2` to
    test this routine.

    Add print statements to the subroutine so that it gives results similar
    to::

        $ make test -f Makefile2
        Wrote data to input_data.txt
        ./test.exe
         n =           20
        Solving tridiagonal system with n =     10
        Solving tridiagonal system with n =      9
        Solving tridiagonal system with n =     10
        Solving tridiagonal system with n =      9
        Computed G0 =   0.1186D+02  G1 =   0.1167D+02  z =  -0.6111D+02
         error_max =   3.99812551484757250E-003


        $ make test -f Makefile2 n=10000
        Wrote data to input_data.txt
        ./test.exe
         n =        10000
        Solving tridiagonal system with n =   5000
        Solving tridiagonal system with n =   4999
        Solving tridiagonal system with n =   5000
        Solving tridiagonal system with n =   4999
        Computed G0 =   0.2442D-01  G1 =   0.2402D-01  z =  -0.6004D+02
         error_max =   1.73004366388340713E-008

    Printing out `G0, G1`, and `z` might be useful for debugging to compare
    with what the Python code produces.

    **NOTE:** The split version gets slightly different results from the direct 
    version, particularly for `n = 10000`.   This tridiagonal matrix has a 
    condition number that grows like :math:`{\cal O}(n^2)`, and so for `n =
    1e5`, perturbations in the way you do the linear algebra can make
    be expected to make relative errors on the order of `1e10` times the 
    machine roundoff, which is comparable to the error in the BVP solution that
    is being computed. (This also suggests that making `n` even larger may
    not give you a better approximation of the true solution, when computing
    in finite precision arithmetic.)

#.  Add another subroutine `solve_bvp_split_omp` to the `bvp_solvers` module
    that has the same calling sequence and produces the same result, but  
    that uses OpenMP in such a way that four different threads make the four
    calls to `solve_bvp_direct`.  

    Do this by using `omp parallel sections
    <https://computing.llnl.gov/tutorials/openMP/#SECTIONS>`_, see for example
    `$UWHPSC/codes/openmp/demo2.f90` or
    `$UWHPSC/codes/adaptive_quadrature/openmp2/adapquad_mod.f90`.

    This will take a bit of thought about what variables should be private
    to each thread and perhaps some rearrangement of the code to make
    sure each thread is solving the desired problem and all four results can
    be combined as needed.  To help debug, you might want to print out
    various things from the serial version of the code and compare to the
    parallel version, and try running with small values of `n`.

    You can call `omp_set_num_threads(4)` in the subroutine and do not
    need to test with a different number of threads.

    **Note:** This is not a great problem for OpenMP since solving a
    tridiagonal system is so quick, and the overhead of forking threads
    will probably make the OpenMP version run slower than the serial version
    unless `n` were very large, but the point is to understand and debug the
    code.

    Provide a new main program `main3.f90` and `Makefile3` that compiles
    with OpenMP and links with OpenMP and the LAPACK libraries, e.g. set:: 

        LFLAGS = -lblas -llapack -fopenmp
        FFLAGS = -fopenmp

    Add print statements to your subroutine so that it gives output such as::

        $ make test -f Makefile3
        test.exe
        Wrote data to input_data.txt
        ./test.exe
         n =           20
         nthreads =            4
        Thread   0 taking from   0.000 to   0.524 with u_mid =   0.000
        Solving tridiagonal system with n =     10
        Thread   1 taking from   0.524 to   1.000 with u_mid =   0.000
        Solving tridiagonal system with n =      9
        Thread   2 taking from   0.000 to   0.524 with u_mid =   1.000
        Solving tridiagonal system with n =     10
        Thread   3 taking from   0.524 to   1.000 with u_mid =   1.000
        Solving tridiagonal system with n =      9
        Computed G0 =   0.1186D+02  G1 =   0.1167D+02  z =  -0.6111D+02
         error_max =   3.99812551484757250E-003

        $ make test -f Makefile3 n=10000
        Wrote data to input_data.txt
        ./test.exe
         n =        10000
         nthreads =            4
        Thread   1 taking from   0.000 to   0.500 with u_mid =   0.000
        Solving tridiagonal system with n =   5000
        Thread   0 taking from   0.500 to   1.000 with u_mid =   0.000
        Solving tridiagonal system with n =   4999
        Thread   2 taking from   0.000 to   0.500 with u_mid =   1.000
        Solving tridiagonal system with n =   5000
        Thread   3 taking from   0.500 to   1.000 with u_mid =   1.000
        Solving tridiagonal system with n =   4999
        Computed G0 =   0.2442D-01  G1 =   0.2402D-01  z =  -0.6004D+02
         error_max =   1.73004366388340713E-008


**Part 2**

In :ref:`lab19` the heat equation will be discussed along with an IPython
notebook illustrating how solutions behave and two numerical methods for
approximating the solution.

For simplicity, we are only considering a special case of the
one-dimensional heat equation :math:`u_t(x,t) = u_{xx}(x,t)` in
which the problem is solved on the interval :math:`0 < x < \pi`,
the boundary conditions are :math:`u(0,t) = u(\pi,t) = 0` for all
:math:`t`, and the initial data is a sine wave of the form `u_0(x)
= \sin(kx)` for some integer `k`.  You might want to experiment
with initial data that is a linear combination of different "Fourier
modes", as illustrated in the notebook.

The directory `$UWHPSC/homeworks/project/part2` contains some files that
implement the explicit method discussed in class.  You can do, for example::

    $ make test -f Makefile1 

and you can vary `n, k, tfinal,` and `nsteps` by specifying at the command
line, e.g. ::

    $ make test -f Makefile1 n=100 k=5 nsteps=500

The main program prints out the max-norm error at the final time and also
produces a file `solution.txt` that contains the approximate and true solution at
the final time.

#. Add a second subroutine to the file `heat_solvers.f90` that implements the
   implicit Crank-Nicolson method that will be discussed in :ref:`lab19`.
   Name this subroutine `solve_heat_implicit` and it should have the same calling
   sequence as the `solve_heat_explicit`.  

   Add a parameter `method` to `main1.f90` so that if `method==1` then the 
   explicit method is used and if `method==2` then the implicit method is used.
   Add this also to `Makefile1` so that a value is written to `input_data.txt`
   and then read by the main program, similar to the other parameters.
   (You can give it the default value 1).
   Note that the two methods do not give the same approximate solution (or
   error), but test that both give results that agree with the IPython notebook.

   To implement this method, you will have to solve a tridiagonal system of
   equations every time step.  You can use the LAPACK routine `dgtsv` (or `dptsv`
   if you prefer).  Note that either of these routines overwrites the input
   arrays that describe the matrix with the LU factorization, so be careful if 
   you are using this in a loop where you have more than one system to solve!

#. Create a new main program `main2.f90` based on your modified `main1.f90` that
   outputs the solution at every time step to a file `frames.txt`.  Use the same
   format as currently used to write to `solution.txt`, but add to the file every
   time step, and also write the initial data before starting to solve the
   problem.  So after running the code the file `frames.txt` should have
   `(nsteps+1)*(n+2)` lines (since each `u` solution vector has `n+2` elements).

   Do not modify the subroutines in `heat_solvers.f90` to do this.  Instead,
   have a loop in the main program that calls `solve_heat_explicit` or
   `solve_heat_implicit` repeatedly, `nsteps` times, taking a single time step
   with each call and then writing the solution before the next call.

#. Write a Python script `animate.py`
   that reads `n` and `nsteps` from `input_data.txt` and 
   reads all the solutions from `frames.txt` and produces an animation in a file 
   `heat.html`.  Use `JSAnimation` and the `JSAnimation_frametools.py` module
   from :ref:`lab15`. 

   In :ref:`lab20` we will look at an example of doing this for a different
   problem, so if you're not sure how to do it, take a look at that Lab.

   Create a `Makefile2` with a phony target `movie` so that you can do, for example,

        $ make movie -f Makefile2  k=4 n=50 nsteps=40 method=1

        $ make movie -f Makefile2  k=4 n=50 nsteps=40 method=2

   and create the animations shown at 

   * `<http://faculty.washington.edu/rjl/heat_explicit.html>`_
   * `<http://faculty.washington.edu/rjl/heat_implicit.html>`_

   illustrating that the implicit method is more stable.



             
To submit
---------

* At the end, you should have committed the following 
  files to your repository:

  **Part 1**

  * `$MYHPSC/project/part1/Makefile1`  (unchanged from original)
  * `$MYHPSC/project/part1/main1.f90`  (unchanged from original)
  * `$MYHPSC/project/part1/problem.f90`  (unchanged from original)
  * `$MYHPSC/project/part1/bvp_solvers.f90` (with 3 subroutines)
  * `$MYHPSC/project/part1/Makefile2`  
  * `$MYHPSC/project/part1/main2.f90` 
  * `$MYHPSC/project/part1/Makefile3`  
  * `$MYHPSC/project/part1/main3.f90` 

  **Part 2**

  * `$MYHPSC/project/part2/problem.f90`  (unchanged from original)
  * `$MYHPSC/project/part2/Makefile1`  
  * `$MYHPSC/project/part2/main1.f90`  
  * `$MYHPSC/project/part2/heat_solvers.f90`  
  * `$MYHPSC/project/part2/Makefile2`  
  * `$MYHPSC/project/part2/main2.f90`  
  * `$MYHPSC/project/part2/animate.py`  (Script to create the animation) 
  * `$MYHPSC/project/part2/JSAnimation_frametools.py`  (optional, unchanged from
    original)


* **Please be sure you have the specified directory and file names.**
  It is hard to grade otherwise, and points will be deducted.
  

* Make sure you push to bitbucket after committing.

* Submit the commit number that you want graded by following the link
  provided on the `Canvas page for the project
  <https://canvas.uw.edu/courses/893991/assignments/2520179>`_.

* Also take this `survey
  <https://canvas.uw.edu/courses/893991/quizzes/782895>`_
  worth 10 points.
  
  Before doing this survey, please first do the course evaluation. This
  is being done on-line this year at the link
  `<https://uw.iasystem.org/survey/128096>`_.
  Anonymous results will be available to the instructor and TAs well after
  the quarter has ended.

* On-campus students:
  UW Libraries and UW Information Technology want to know what students
  think about the Active Learning Classrooms in Odegaard Library. Your
  feedback will shape recommendations for new classrooms. Take their brief
  survey here: `<https://catalyst.uw.edu/webq/survey/fournier/239147>`_.







