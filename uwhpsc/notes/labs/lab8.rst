

.. _lab8:

Lab 8: Thursday April 24, 2014
=============================

The code for this lab can be found in `$UWHPSC/labs/lab8`.

Lab started with a demo of using compiler flags and `gdb` to debug Fortran
code, using the code in  `$UWHPSC/labs/lab8/demo`.

Running it without any compiler flags gives no error by there is a `NaN`
value in the results::

    $ cd $UWHPSC/labs/lab8
    $ make run

     The max value of y is   0.99973335466585400     
     x(501) is    0.0000000000000000     
     y(501) is                        NaN

Running it with compiler flags to catch floating point exceptions::

    $ make debug
    make: *** [debug] Floating point exception: 8

Once it's compiled with the flags specified in the Makefile for the `debug`
target, the debugger `gdb` can be used to run
the code and figure out where it died::

    $ gdb ./a.out
    (gdb) run

    Program received signal SIGFPE, Arithmetic exception.
    0x0000000000400a6d in demo () at demo.f90:12
    12          y(j) = sin(x(j)) / x(j)

    (gdb) p j
    $1 = 501

    (gdb) p x(j)
    $2 = 0

Many commands are available in `gdb`, see for example 
`this documentation <http://www.sourceware.org/gdb/current/onlinedocs/gdb.html>`_.

**Note:** On the Mac with the Mavericks OS, `gdb` has been replaced by
`lldb`.  See `<http://lldb.llvm.org/index.html>`_ for more information.

The page `GDB TO LLDB COMMAND MAP <http://lldb.llvm.org/lldb-gdb.html>`_
gives a good summary of both `gdb` and `lldb` commands and the relation
between them.


Debugging compile-time errors
-----------------------------

The code in `$UWHPSC/labs/lab8/problem1` does not compile.  See if you can
find and fix all the errors.  See `$UWHPSC/labs/lab8/problem1b`
for a corrected version (and see the `README.txt` file in that directory for
some comments).

Debugging run-time errors
-----------------------------

The code `$UWHPSC/labs/lab8/problem2/array1.f90` does not run properly.  
See if you can find and fix all the errors.  See
`$UWHPSC/labs/lab8/problem2/array1b.f90` for a corrected version and use
`diff` to see the differences.
