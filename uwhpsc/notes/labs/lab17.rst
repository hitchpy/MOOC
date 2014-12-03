

.. _lab17:

Lab 17: Tuesday May 27, 2014
=============================

Announcements
-------------

* Due tonight: Homework 4, lecture quizzes, lab quiz.

* Office hours today: 
  
 * In Odegaard after class, and then in Lewis 328 until 5pm,
 * Online 5-6pm, look for announcement on Canvas, but may use GoToMeeting

* Final project:  Part 1 will be posted soon, discussed on Thursday.

  Due Wednesday, June 11.  (But don't wait to the last minute!)

Today's lab
----------------

* Go through the notebook `$UWHPSC/labs/lab17/Tridiagonal.ipynb`, also 
  visible at
  `<http://nbviewer.ipython.org/url/faculty.washington.edu/rjl/notebooks/Tridiagonal.ipynb>`_.
   
* Work in pairs on writing a Fortran program to solve the same 
  tridiagonal linear system considered in the notebook with::

      A =
        [[   1.  200.    0.    0.    0.]
         [  10.    2.  300.    0.    0.]
         [   0.   20.    3.  400.    0.]
         [   0.    0.   30.    4.  500.]
         [   0.    0.    0.   40.    5.]]

      b = 
        [ 201.  312.  423.  534.   45.]

* Use the LAPACK surbroutine `dgtsv`.  See the 
  documentation at `<http://www.netlib.no/netlib/lapack/double/dgtsv.f>`_.

* **There is quiz for Lab 17**
