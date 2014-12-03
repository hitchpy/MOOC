
.. _homework2:

==========================================
Homework 2 [2014]
==========================================


.. warning:: Assignment modified in Part 9.

Due Thursday, April 24, 2014, by 11:00pm PDT.
See below for how to submit.


The goals of this homework are to:

* Learn more about "best practices in scientific computing".

* Get some more experience with Python and *numpy*.
* Get some experience with Fortran code and Makefiles.


Before tackling this homework, you should read some of the class notes and
links they point to.  In particular, the following sections are relevant:

* :ref:`python`
* :ref:`numerical_python`
* :ref:`special_functions`
* :ref:`fortran`
* :ref:`makefiles`

#.  **Reading assignment**

    Download and read the paper `Best Practices for Scientific Computing 
    <http://www.plosbiology.org/article/info%3Adoi%2F10.1371%2Fjournal.pbio.1001745>`_ by
    G. Wilson, D. A. Aruliah, C. T. Brown, et. al.
    There is a link to download the pdf file on the page linked, and the
    article is open access so anyone should be able to download it.

    You should be pleased to find that you are now starting to follow many
    of these best practices, but there are many good tips in the paper that
    have not been covered in lectures.

    There will be a quiz on this paper to complete as part of the homework
    assignment.  This quiz can be found on the Canvas page
    `Homework 2 quiz <https://canvas.uw.edu/courses/893991/quizzes/779589>`_.

#. **Survey assignment**

   Please take `this survey <https://canvas.uw.edu/courses/893991/quizzes/779591/>`_
   (worth 5 points) to give us some feedback on how the class is going
   so far.
    
#.  **Programming assignments**

    You should create a new subdirectory `homework2` (of the same private
    repository you have used for submitting previous  homeworks).  
    Develop your code for the problems below in this directory
    and feel free to commit as
    often as you like, it will help you recover from blunders.

    Note that to grade this homework, we will try using the functions
    you write with different input values than the ones used in the examples 
    below, so you might want to test your functions for other reasonable
    input values.

#.  Recall that computer hardware can only do basic arithmetic, so if you
    want to evaluate some special function such as square root or a
    trigonometric function, some algorithm must be implemented by someone to
    compute an approximation for arbitrary input values.  To approximate
    :math:`\sqrt{y}`, we saw that one approach is to use Newton's method 
    to approximate a zero of the function :math:`x^2 - y = 0`, for example.
    (See also :ref:`special_functions` for more about this.)

    Another approach is to use a Taylor series approximation, see e.g. the
    `Wikipedia Taylor series page <http://en.wikipedia.org/wiki/Taylor_series>`_
    for a review.  

    Recall that the Taylor series for the exponential function (expanded
    about :math:`x_0 = 0`) is given by

        :math:`\exp(x) = \sum_{j=0}^\infty \frac{x^j}{j!} 
        = 1 + x + \frac 1 2 x ^2 + \frac 1 6 x^3 + \cdots`

    with the convention that :math:`x^0 = 0!  = 1`.


    Create a Python file `taylor.py` in the `homework2` directory of your
    repository and write a function `exp1(x,n)` in this file that 
    returns an approximation to :math:`\exp(x)` based on the Taylor series
    of degree `n`.  (Note that this degree `n` polynomial is obtained by
    truncating the infinite series after `n+1` terms.)


    If you have done this properly, then 
    if you start the Python or Ipython shell in the same directory,
    you should be able to do, for example::

        >>> import taylor
        >>> taylor.exp1(1.,0)
        1.0
        >>> taylor.exp1(1.,1)
        2.0
        >>> taylor.exp1(1.,2)
        2.5
        >>> taylor.exp1(1.,20)
        2.7182818284590455

    

    **Hint:** Note that term of degree `j` in the series can be computed 
    from the previous term by multiplying by `x` and dividing by `j`. If you
    use this trick you won't need the factorial function `math.factorial`.
    See :ref:`fortran_taylor` and
    `$UWHPSC/codes/fortran/taylor.f90` for an example of this same
    idea used in a Fortran version.

#.  Add some debugging statements to your function, with an optional argument
    `debug` with the default value `False` (so the examples above still give
    the same output) but so that setting `debug=True` causes output similar to
    this::

        >>> taylor.exp1(1.,5,debug=True)
        j = 1, term = 1.000000000000000e+00
          partial_sum updated from 1.000000000000000e+00 to 2.000000000000000e+00
        j = 2, term = 5.000000000000000e-01
          partial_sum updated from 2.000000000000000e+00 to 2.500000000000000e+00
        j = 3, term = 1.666666666666667e-01
          partial_sum updated from 2.500000000000000e+00 to 2.666666666666667e+00
        j = 4, term = 4.166666666666666e-02
          partial_sum updated from 2.666666666666667e+00 to 2.708333333333333e+00
        j = 5, term = 8.333333333333333e-03
          partial_sum updated from 2.708333333333333e+00 to 2.716666666666666e+00
        2.7166666666666663

    You probably won't need these statements for this function,  but similar
    statements might be useful in the next part.


#.  Create a Python function `sin1(x,n)`  (in the same file
    `taylor.py` as the function `exp`) that approximates the sine
    function at a point `x` by evaluating the Taylor series approximation of
    degree `n`.  Use the Taylor series expansion about :math:`x_0=0`,
    also known as the Maclaurin series:

        :math:`\sin(x) = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \frac{x^7}{7!}
        + \cdots`

    Note that the degree 5 and 6 approximations only have three nonzero terms, 
    the degree 7 and 8 approximations have four nonzero terms, etc.  

    You should get results like::

        >>> taylor.sin1(pi/2, 2)
        1.5707963267948966

        >>> taylor.sin1(pi/2, 3)
        0.9248322292886504

        >>> taylor.sin1(pi/2, 4)
        0.9248322292886504

        >>> taylor.sin1(pi/2, 5)
        1.0045248555348174

    Add a debug option as in `exp1`.  

    **Hint:**  You might find it convenient to have a variable `term` that
    is updated as for the exponential function but then is multiplied by
    `s` before adding in to the partial sum, where `s` takes the appropriate
    value :math:`+1,~-1,` or 0 depending on `j`.

#.  See what happens if you call your function `exp1` or `sin1` with
    negative values of `n`, or with non-integer real numbers.  Add
    some input-checking to each function so that a non-negative integer
    value of `n` is required.  If an invalid value is detected, print an
    error message and return the special value `numpy.nan` ("not a number",
    similar to the Matlab `NaN`).  For example::

        >>> taylor.exp1(1., -3)
        *** Invalid input -- n must be non-negative integer
        nan

#.  The code `$UWHPSC/codes/fortran/taylor.f90` contains a main program and
    subroutine for approximating exp(x) by a Taylor series.  Split this code
    up into two separate files `taylor_main.f90` and `exptaylor.f90` and add
    a Makefile based on `$UWHPSC/codes/fortran/multifile1/Makefile5`
    so that you can do::

        $ make exp_output.txt

    to create a file `exp_output.txt` containing::

         x =    1.0000000000000000     
         n =           20
         exp_true  =    2.7182818284590451     
         exptaylor =    2.7182818284590455     
         error     =   4.44089209850062616E-016

    **Note:** Also modify the main program so it prints the value of `n`
    as shown above.



#.  **Only 583 students need to do this part**

    (483 students are encouraged to do these parts too, 
    but they will not count towards the score -- the parts will be weighted
    differently for 583 students.
    Note that undergrads who registered for 583A will be treated as 483
    students.)

    Add another file `sinetaylor.f90` that computes the approximation to the 
    sine function, as you did in the Python version.  
    
    .. warning:: Assignment modified here:

    Also create a new main program `taylor_main2.f90`
    that calls this subroutine.  Add these to the *same* Makefile so that ::

        $ make sine_output.txt

    gives something sensible and `make exp_output.txt` still gives the
    previous results.


#.  **Only 583 students need to do this part**

    The `Wikipedia Taylor series page <http://en.wikipedia.org/wiki/Taylor_series>`_
    shows a nice plot of Taylor series
    approximations to the sine function for different orders.  The gnuplot
    commands that created this plot can be found at
    `http://en.wikipedia.org/wiki/File:Sintay_SVG.svg
    <http://en.wikipedia.org/wiki/File:Sintay_SVG.svg>`_.

    Write a script `plot_taylor.py` to produce a *similar* plot that shows
    the sine function and approximations for `n = 1,3,5,7` over the same range
    of `x` values.  You don't need to try to match the colors or add the
    grid lines.

To submit
---------

* At the end, you should have committed the following 
  files to your repository:

  * `$MYHPSC/homework2/taylor.py`
  * `$MYHPSC/homework2/taylor_main.f90`
  * `$MYHPSC/homework2/exptaylor.f90`
  * `$MYHPSC/homework2/Makefile`

  583 students should also have the files

  * `$MYHPSC/homework2/taylor_main2.f90`
  * `$MYHPSC/homework2/sinetaylor.f90`
  * `$MYHPSC/homework2/plot_taylor.py`
    
    You do not need to submit the png file of the figure this creates.


  Make sure you push to bitbucket after committing.

* Submit the commit number that you want graded by following the link
  provided on the `Canvas page for Homework 2
  <https://canvas.uw.edu/courses/893991/wiki/homework-2>`_.
  If you submit the wrong thing or make further changes to your work
  before the due date, you can simply resubmit new information at the same
  link.

