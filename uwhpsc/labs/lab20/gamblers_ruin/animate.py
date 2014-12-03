import matplotlib
matplotlib.use('Agg')

from pylab import *
import JSAnimation_frametools as J

plotdir = '_plots'  # to store png files for each figure
J.make_plotdir(plotdir, clobber=True)  # ok to clobber if it already exists


input_data = loadtxt('input_data.txt')
p = float(input_data[2])
seed = int(input_data[4])

n1, n2 = loadtxt('walk.txt', unpack=True)

nsteps = len(n1)
nmax = max(n1.max(), n2.max())  # to set axis limits


for k in range(nsteps):
    steps = range(0,k+1)
    n1steps = n1[0:k+1]
    n2steps = n2[0:k+1]
    clf()
    plot(steps,n1steps,'bo-',label='n1')
    plot(steps,n2steps,'ro-',label='n2')
    axis([0,nsteps,0,nmax])
    legend(loc = 'lower left')
    title('seed = %s' % seed)

    J.save_frame(k, plotdir, verbose=True)

anim = J.make_anim(plotdir)

description = "Gambler's ruin random walk with p = %5.3f" % p

J.make_html(anim, file_name="RandomWalk.html", title=description)
    
