PyGMIN, a python version of GMIN, OPTIM, and PATHSAMPLE.  GMIN was originally
written by David Wales and is designed to find the global energy minimum for
molecular systems.

The code is hosted at
https://github.com/js850/PyGMIN

for documentation, see
http://js850.github.com/PyGMIN/

The code project upon which this python package is based can be found at
http://www-wales.ch.cam.ac.uk/software.html


INSTALLATION
============

Required packages
-----------------

for compilation::
  cython
  f2py (comes bundled with numpy, see below)

python packages::

  numpy (with f2py)
  scipy
  networkx
     https://networkx.lanl.gov
     for graph functionality
  SQLAlchemy 0.7 
     http://www.sqlalchemy.org/
     for managing database of stationary points 


  hungarian (optional, use pip/easy_install to get newest version with bugfixes!)
     for min dist routines
  scikits.sparse (optional, for use of sparse Cholesky decomposition methods when
     calculating rates)
     
All the above packages can be installed via the python package managers pip or
easy_install.  However, some of the packages (numpy, scipy) have additional
dependencies and it can be more convenient to use the linux package manager
(apt, yum, ...).

If you want to use the gui you need:
  matplotlib 
  qt4 python bindings
  opengl python bindings

The Ubuntu packages are (apt-get)::

  python-qt4
  python-opengl 
  python-qt4-gl

  pymol (optional)

Or Fedora (yum)::
   PyQt4
   PyOpenGl


installing prerequisites on Ubuntu
----------------------
if you're running ubuntu, you can get all the prerequisites with the following
  commands

$ sudo apt-get install python-numpy python-scipy python-matplotlib python-qt4 python-opengl python-qt4-gl python-pip cython pymol
$ pip install --user networkx sqlalchemy hungarian
(in the above, the flag --user will install localy.)


Compilation
-----------

Compilation is required for use of the fast potentials, those written in C++
and/or fortran).  Theoretically you should be able to use any fortran compiler,
but we mostly use gfortran, so it's the least likely to have problems.  The
package uses the standard python setup utility (distutils).  There are lots of
options for how and where to install. For more information::
  
  python setup.py --help 
  python setup.py --help-commands

Developers probably want to install "in-place", i.e. build the extension
modules in their current directories::

  python setup.py build_ext -i --fcompiler=gfortran

Users can install pygmin in the standard python package location::

  python setup.py build --fcompiler=gfortran
  python setup.py install [--user]

where --user installs it in $HOME/.local/


PYTHONPATH  
----------
make sure to add the install directory to your PYTHONPATH environment variable.
This is not necessary if you install to a standard location.


Installing on Mac
-----------------
When I used macports, everything installed very easily.  However, at least on
my machine the OpenGL / Qt implementation was broken and I could never get it
to work.  So I aborted using macports and installed everything by hand.  I
eventually got (almost) everything working.  Here are a few hints

use the Enthought python distribution instead of the prepackaged one.  This
seems to be standard, plus it includes numpy and scipy
http://www.enthought.com/products/epd.php

if you want to use the gui you have to install pyqt and its dependencies.  This
is not as simple as it should be.  I roughly followed the instructions here
http://www.noktec.be/python/how-to-install-pyqt4-on-osx .  Even though my mac is
64 bit I had to compile things with --arch=i386.  I think I even had to install
Qt from source to get it with the 32 bit architecture.    Also, when compiling
pyqt, you have to point to the qmake in the Qt install folder, not
/usr/bin/qmake as the website says.  Sorry for the vague instructions, I'm
writing this from memory.  If you have more complete installation help please
email or submit a pull request.

Installing pymol was easy with macports, but since I've abandoned macports 
I haven't gotten it installed and running (I haven't tried the pay version) 


=== Running ===
The documentation can be found at

http://js850.github.com/PyGMIN/

Alternatively the documentation can be compiled using sphinx ( in doc/ run $ make html ).

Also, see the examples in the examples folder.
