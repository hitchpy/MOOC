
.. _homework4_solution:

==========================================
Homework 4 [2014] Discussion of solution
==========================================

Sample solutions can be found in `$UWHPSC/solutions/homework4`.

**Reminder to help debugging Fortran codes:**

*   Use the compiler flag `-g` so that you can run the executable
    with GDB (or LLDB on Mac
    Mavericks), which might help determine where the code dies if
    it does, and can help in examining variables as it runs by
    setting breakpoints.

*   Additional compiler flags are very helpful, e.g. `-fbounds-check`
    to have it check
    array bounds, so that for example if `x` is dimensioned as
    `x(0:10)` and you try to access `x(11)` then it will give an
    error rather than taking whatever garbage is in the element of
    memory that this points to.

*   See :ref:`fortran_debugging` and :ref:`gfortran_flags` for more tips.

*   With the Makefiles we are using, you can set e.g.::

        FFLAGS = -g -Wall -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid

    or even set these on the command line with something like::

        $ make test1 -f Makefile1 FFLAGS="-g -Wall -fbounds-check"

*   Remember that to have these take effect you will have to remove
    all the old `.o` files
    first so that it recompiles the code rather than thinking they
    are up to date.  You can do `rm *.o` or, if the Makefile is set
    up with a `clean` target, something like::

        $ make clean -f Makefile1



**Some common errors that people made:**

**Part 1**

* The most common error was blindly following the coarse grain
  example and dividing up `n` instead of `n_samples` between threads.

* Incorrect OpenMP reductions.

**Part 2**

* A variety of problems coding `main2f.90`.

