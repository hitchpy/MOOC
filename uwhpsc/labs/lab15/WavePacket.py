
import matplotlib
matplotlib.use('Agg') # graphics backend for hardcopy
from pylab import *

from JSAnimation import IPython_display
import JSAnimation_frametools as J

plotdir = '_plots'  # to store png files for each figure
J.make_plotdir(plotdir, clobber=True)  # ok to clobber if it already exists

x = linspace(0., 10., 1001)
nsteps = 21
dt = 11. / (nsteps+1)

plotdir = '_plots'

for n in range(nsteps+1):
    t = n*dt
    x0 = x - 1. - t
    u = sin(20.*x0) * exp(-5.*x0**2)
    clf()
    plot(x,u,'b')
    ylim(-1.2, 1.2)
    title('Wave packet at t = %8.4f' % t, fontsize=20)
    
    # Save this frame:
    J.save_frame(n, plotdir,verbose=False)
    #show()

anim = J.make_anim(plotdir)

J.make_html(anim, file_name="Wave.html", title="Advecting wave packet")
