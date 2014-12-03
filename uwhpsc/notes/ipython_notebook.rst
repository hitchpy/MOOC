

.. _ipython_notebook:

=============================================================
IPython_notebook
=============================================================

The IPython notebook is fairly new and changing rapidly.  The version
originally installed in the class VM is version 0.10.  To get the latest
development version, which has some nicer features, do the following::

    $ cd
    $ git clone https://github.com/ipython/ipython.git
    $ cd ipython
    $ sudo python setup.py install

Then start the notebook via::

    $ ipython notebook --pylab inline

in order to have the plots appear inline.  If you leave off this argument, a
new window will be opened for each plot.

Read more about the notebook in the `documentation
<http://ipython.org/ipython-doc/dev/interactive/htmlnotebook.html>`_

See some cool examples in the `IPython notebook viewer
<http://nbviewer.ipython.org/>`_.

See also :ref:`sagemath`.

Interactive notebooks
---------------------

`IPython 2.0 <http://ipython.org/ipython-doc/dev/whatsnew/version2.0.html>`_
(released April 1, 2014) includes `interactive widgets <http://nbviewer.ipython.org/github/ipython/ipython/blob/2.x/examples/Interactive%20Widgets/Index.ipynb>`_

See these `Tips for installing IPython 2.0 <http://ipython.org/install.html>`_
on your own computer.

SageMathCloud does not yet have IPython 2.0 and for various technical
reasons will not have it for a while.

SageMathCloud does now have `mpld3 <https://github.com/jakevdp/mpld3>`_, which
allows zooming in on plots.  For a demo, see :ref:`lab13` or 
`Jake's blog post on mpld3
<http://jakevdp.github.io/blog/2013/12/19/a-d3-viewer-for-matplotlib/>`_

In addition to `mpld3`, Jake Vanderplass has also developed:

* `ipywidgets <https://github.com/jakevdp/ipywidgets>`_ similar to the 2.0
  widgets in some ways but persistent also if the "static" notebook is
  viewed via `nbviewer <http://nbviewer.ipython.org>`_.
  See `Jake's blog post on ipywidgets
  <http://jakevdp.github.io/blog/2013/12/05/static-interactive-widgets/>`_.

* `JSAnimation <https://github.com/jakevdp/JSAnimation>`_ for persistent
  animations.  This is used for example for all the animations in the
  `Clawpack galleries <http://clawpack.github.io/doc/galleries.html>`_.
  (`Clawpack <http://www.clawpack.org>`_ is an open source software 
  package developed by Randy LeVeque and many others for solving 
  hyperbolic partial differential equations.)

