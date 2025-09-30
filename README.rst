..
    .. image:: https://travis-ci.org/farrelljd/pele.svg?branch=modernise
        :target: https://travis-ci.org/pele-python/pele?branch=master

    .. image:: https://coveralls.io/repos/pele-python/pele/badge.png?branch=master
        :target: https://coveralls.io/r/pele-python/pele?branch=master

    .. image:: https://landscape.io/github/pele-python/pele/master/landscape.svg
       :target: https://landscape.io/github/pele-python/pele/master
       :alt: Code Health

pele : Python Energy Landscape Explorer
+++++++++++++++++++++++++++++++++++++++

Tools for global optimization and energy landscape exploration.

Source code: https://github.com/pele-python/pele

Documentation: http://pele-python.github.io/pele/

.. figure:: ./images/lj38_gmin_dgraph.png

  Images: The global minimum energy structure of a 38 atom Lennard-Jones cluster.  On
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

Installation
============

Install via conda and pip. In this directory:

.. code-block:: bash

    conda env create -f environment.yml
    conda activate pele
    pip install --no-build-isolation .

Works like a dream on Ubuntu 22.04---let us know if it works on your system.

Running
=======

You can find examples of how to run pele in the examples folder.  More
information can be found in the documentation at

http://pele-python.github.com/pele/


Notes
=====
In the long-long-ago, pele was known as pygmin.

Tests
=====
Pele has a large suite of unit tests.  They can be run using the nose testing
framework (which can be installed using pip).  The tests are run from the top
directory with this command::

  nosetests pele

We also have test for our c++ code writen in pure c++.  These are stored in
the folder cpp_tests/ and can be compiled using CMake.
