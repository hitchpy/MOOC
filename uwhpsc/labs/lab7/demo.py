"""
Sample code to demo debugging in IPython with pdb
"""

import numpy as np
# np.seterr(all='raise')  # raise floating point exceptions

def printmax():
    x = np.linspace(-20, 20, 1001)
    y = np.zeros(x.shape)
    for j in range(len(x)):
        y[j] = np.sin(x[j])/x[j]

    print "the max value of y is ", y.max()
    return y

if __name__=="__main__":
    y = printmax()

