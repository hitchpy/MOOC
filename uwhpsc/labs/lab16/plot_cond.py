
import matplotlib
matplotlib.use('Agg')  # backend
from matplotlib import pyplot as plt
import numpy as np

try:
    data = np.loadtxt('cond.txt')
except:
    print "*** Error, could not read file cond.txt with loadtxt"
    data = None

if data is not None:
    n = data[:,0]     # first column
    cond = data[:,1]  # second column
else:
    n = np.linspace(1,20,20)

cond_formula = 0.01133 * np.exp(3.49*n)

plt.figure()

# semilog plot:
plt.semilogy(n,cond_formula,'ko-',label='from formula')
if data is not None:
    plt.semilogy(n,cond,'rx-',label='estimated')
plt.legend(loc='upper left')
plt.title("Condition number of n by n Hilbert matrix")
plt.xlabel('n')

fname = 'cond.png'
plt.savefig(fname)
print "Created ",fname



