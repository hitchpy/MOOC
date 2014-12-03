
.. _animation:

Animation in Python
===================

`matplotlib` has a package `animation` that can be used directly, 
see for example
`<http://matplotlib.org/examples/animation/dynamic_image2.html>`_
or `this blog post
<http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/>`_.

Nicer webpage animations (within IPython notebooks or as stand-alone movies)
can be created using the package `JSAnimation` created by
Jake Vanderplas.  For an example see this `rendered example
<http://nbviewer.ipython.org/github/jakevdp/JSAnimation/blob/master/animation_example.ipynb>`_
or :ref:`lab15`.

Installing JSAnimation
======================

First clone it from Github::

    $ cd $HOME
    $ git clone https://github.com/jakevdp/JSAnimation.git
    $ cd JSAnimation
    
On your own laptop or the VM, you can probably install it via::
   
    $ python setup.py install

On SageMathCloud you don't have access to the system folder where it
normally installs Python packages, but you can install it for use in
one project only via::
   
    $ sage -python setup.py install --user

Then you should be able to open Python or IPython and `import JSAnimation`
       
Alternative to installing
-------------------------

Rather than installing it as a package, 
you can just add the JSAnimation directory to 
the Python search path.  For running it from scripts in a bash shell, 
you might want to add this line to your 
Or you can just add $HOME/JSAnimation to your PYTHONPATH environment
variable::

    $ export PYTHONPATH=$PYTHONPATH:$HOME/JSAnimation

This appends the path to the end of whatever path is already specified in
this environment variable.

As a last resort, you can also modify the path from within a Python
session::

    >>> import os, sys
    >>> HOME = os.environ['HOME']
    >>> JSAnimation_path = os.path.join(HOME, 'JSAnimation')
    >>> sys.path.append(JSAnimation_path)


JSAnimation_frametools
----------------------

For animations of complex plots, it is sometimes easier to simply
create the plot for each frame as usual using `matplotlib` and then
save a `.png` file for each frame.  These can be created and then
combined to create an animation using some tools in
`JSAnimation_frametools.py`, currently found in :ref:`lab15` along
with some demos.

Matplotlib issues
------------------

JSAnimation requires a recent version of `matplotlib`.
(In particular, older Ubuntu versions may not have a recent version.)
If you're having problems with `matplotlib` in this context, you might
want to try using the `Anaconda Python distribution
<https://store.continuum.io/cshop/anaconda>`_, or switch to :ref:`smc`.

