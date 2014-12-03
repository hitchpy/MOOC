import matplotlib
matplotlib.use('Agg')

from pylab import *
import JSAnimation_frametools as J

plotdir = '_plots'  # to store png files for each figure
J.make_plotdir(plotdir, clobber=True)  # ok to clobber if it already exists

input_data = loadtxt('input_data.txt')
n = int(input_data[0])
nterms = int(input_data[1])

x, f = loadtxt('xf.txt', unpack=True)

data = loadtxt('frames.txt')

terms = reshape(data[:,0], (n+2, nterms), order='F')
fsums = reshape(data[:,1], (n+2, nterms), order='F')

for k in range(nterms):
    term = terms[:,k]
    fsum = fsums[:,k]

    clf()
    subplot(2,1,1)
    plot(x,fsum,'b-')
    plot(x,f,'r-')
    title('Sum of first %s terms in Fourier series' % (k+1))
    axis([0,pi,-1.0,30.0])

    subplot(2,1,2)
    plot(x,term,'b-')
    title("k = %s term in Fourier series" % (k+1))
    plot(x, 0*x, 'k-')  # plot x-axis
    axis([0,pi,-8,10])


    J.save_frame(k, plotdir, verbose=True)

anim = J.make_anim(plotdir)

description = "Fourier sine series for the function exp(x)"

J.make_html(anim, file_name="FourierSum.html", title=description)
    
