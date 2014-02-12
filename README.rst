pele : Python Energy Landscape Explorer
+++++++++++++++++++++++++++++++++++++++

Tools for global optimization and energy landscape exploration.

Source code: https://github.com/pele-python/pele

Documentation: http://pele-python.github.io/pele/



.. figure:: lj38_gmin_dgraph.png

  The global minimum energy structure of a 38 atom Lennard-Jones cluster.  On
  the right is a disconnectivity graph showing a visualization of the energy
  landscape.  The competing low energy basins are shown in color.

pele is a python partial-rewriting of GMIN, OPTIM, and PATHSAMPLE: fortran
programs written by David Wales of Cambridge University and collaborators
(http://www-wales.ch.cam.ac.uk/software.html).  

Description
===========
pele has tools for energy minimization, global optimization, saddle point
(transition state) search, data analysis, visualization and much more.  Some of
the algorithms implemented are:

1. Basinhopping global optimization

#. LBFGS minimization (plus other minimizers)

#. Single ended saddle point search:

  - Hybrid Eigenvector Following

  - Dimer method

4. Double ended saddle point search

  - Nudged Elastic Band (NEB)

  - Doubly Nudged Elastic Band (DNEB)

5. Disconnectivity Graph visualization

6. Structure alignment algorithms

7. Thermodynamics (e.g. heat capacity) via the Harmonic Superposition Approximation

8. Transition rates analysis

INSTALLATION
============

Required packages
-----------------

for compilation:

  1. fortran compiler

  #. c / c++ compiler

python packages:

  1. numpy: 
       We use numpy everywhere for doing numerical work.  It also installs f2py which
       is used to compile fortran code into modules callable by python.

  #. scipy:
       For some of the optimizers and various scientific tools

  #. networkx: 
       For graph functionality. https://networkx.lanl.gov

  #. matplotlib:
       For making plots (e.g. disconnectivity graphs)

  #. SQLAlchemy 0.7: 
       For managing database of stationary points.  http://www.sqlalchemy.org/

  #. hungarian: 
       For permutational alignment

  #. pyro4: 
       For parallel jobs

  #. scikits.sparse: optional 
       for use of sparse Cholesky decomposition methods when calculating rates

  #. pymol: optional
       for viewing molecular structures


All the above packages can be installed via the python package manager pip (or
easy_install).  However, some of the packages (numpy, scipy) have additional
dependencies and it can be more convenient to use the linux package manager
(apt, yum, ...).

If you want to use the gui you will additionally need:

  1. qt4 and qt4 python bindings

  #. opengl python bindings

  The Ubuntu packages (apt-get) for these are: python-qt4, python-opengl, and python-qt4-gl

  In fedora Fedora (yum) you will want the packages: PyQt4, and PyOpenGl


Installing prerequisites on Ubuntu
----------------------------------
if you're running ubuntu, you can get all the prerequisites with the following
commands::

  $ sudo apt-get install python-numpy python-scipy python-matplotlib python-qt4 python-opengl python-qt4-gl python-pip cython pymol
  $ pip install --user networkx sqlalchemy hungarian pyro4 brewer2mpl

(in the above, the flag --user will install localy.)


Compilation
-----------

Compilation is required for use of the fast potentials, those written in C
and/or fortran.  Theoretically you should be able to use any fortran compiler,
but we mostly use gfortran, so it's the least likely to have problems.  This
package uses the standard python setup utility (distutils).  There are lots of
options for how and where to install. For more information::
  
  $ python setup.py --help 
  $ python setup.py --help-commands

Developers probably want to install "in-place", i.e. build the extension
modules in their current directories::

  $ python setup.py build_ext -i --fcompiler=gfortran

Users can install pele in the standard python package location::

  $ python setup.py build --fcompiler=gfortran
  $ python setup.py install [--user]

where --user installs it in $HOME/.local/


PYTHONPATH  
----------
If you do an in-place install, make sure to add the install directory to your
PYTHONPATH environment variable.  This is not necessary if you install to a
standard location.


Installing on Mac
-----------------

Everything installed very easily on my Macbook Air OSX Version 10.75 except the
things needed for the gui.  There is a problem (not related to pele) with the
combination of PyQt4, Qt4, and OpenGL.  If you don't want the gui you should be
golden, but if you do, you may have to install a few things from source.  Below
are the steps I took to get everything working

I use the Enthought python distribution instead of the prepackaged one.  This
seems to be standard, plus it includes numpy and scipy
http://www.enthought.com/products/epd.php

If you want to use the gui you have to install PyQt4 and its dependencies.
This is not as simple as it should be.  Even though my mac is 64 bit I had to
compile everything with --arch=i386.  I even had to install Qt from source to
get it with the 32 bit architecture.   Here are some rough instructions adapted
from http://www.noktec.be/python/how-to-install-pyqt4-on-osx .  That website
gives a good start, but it is not complete.

1. install Qt4.8 from source.  We cannot use the dmg file becuse we need to
   install it for i386 architecture.  
   http://download.qt-project.org/official_releases/qt/4.8/4.8.5/qt-everywhere-opensource-src-4.8.5.tar.gz

   In the directory you unpack the tar.gz file run the following commands.
   http://qt-project.org/doc/qt-4.8/install-x11.html .

   ::

     ./configure -arch i386
     make
     make install

   Make a note of the location of the qmake file that this installs.  We
   will need it for the PyQt4 installation.
  
2. install SIP from source.
   http://www.riverbankcomputing.co.uk/software/sip/download

   In the directory you unpack the tar.gz file run the following commands
   ::

     python configure.py --arch i386
     make
     make install

   This will install SIP for the version of python you use to run configure.py,
   so make sure you're using the correct python version.  Running python
   configure.py --help will tell you which python directory it will be
   installed to.  This should be the same as when you type `which python`
   
3. install PyQt4 from source
   http://www.riverbankcomputing.co.uk/software/pyqt/download .

   In the directory you unpack the tar.gz file run the following commands
   ::

     python configure.py -q <path to qmake in Qt4 folder>  --use-arch i386
     make
     make install

   You must specify (I think) the qmake file that was installed along with Qt4.
   It should be in the Qt4 install directory.

   The same warning for which version of python you use to run configure.py
   applies here as well.

If you have updates or more complete installation instructions please email or
submit a pull request.

Running
=======

You can find examples of how to run pele in the examples folder.  More
information can be found in the documentation at

http://pele-python.github.com/pele/


Notes
=====
pele has recently been renamed from pygmin
