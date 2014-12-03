
"""
Define a function and then print a table of function values.
"""

import numpy as np

def myabs(x):
    """
    This function is written in a complicated way but the formulas
    are meant to reduce to the absolute value function, myabs(x) = |x|.
    """

    y1 = x**2
    y2 = np.sqrt(y1)
    y3 = 1. + np.exp(y2-2.)
    y4 = np.log(y3-1.) + 2.
    y5 = y4 - 3.*(y4 - 2.)
    y6 = (6.-y5)/2.

    return y6

def test_myabs():
    """
    Check a few values against what we expect.
    """
    for x in np.linspace(-2,2,11):
        error_message = "*** Oops, for x = %g, myabs = %g, expected %g" \
            % (x,myabs(x),abs(x)) 
        assert abs(myabs(x) - abs(x)) < 1e-15, error_message


def make_table():
    print "           x                       myabs(x)"
    for x in np.linspace(-2,2,11):
        print "%23.15e   %23.15e  " % (x,myabs(x))


if __name__=="__main__":
    make_table()
    test_myabs()
