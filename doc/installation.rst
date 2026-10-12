.. _installation:

Installation
============

Install from PyPI
-----------------

pele wheels support Python 3.11–3.14 on Linux x86_64 with glibc 2.28 or
newer, and macOS 14 or newer on Intel and Apple Silicon. Windows is not
supported. Install pele into your Python environment with:

.. code-block:: console

   pip install pele

Check that the compiled potential and optimization modules import successfully:

.. code-block:: console

   python -c "from pele.potentials import LJ; from pele.optimize import lbfgs_py; print('pele is ready')"

Then follow :doc:`getting_started` to minimize a configuration and store minima
from a basinhopping run.

Build from source
-----------------

Building from source requires Python 3.11 or newer, C/C++/Fortran compilers,
OpenMP, SUNDIALS in double precision, Eigen, and LAPACK headers. A conda
environment can supply these dependencies:

.. code-block:: console

   conda create -n pele -c conda-forge python=3.12 pip compilers sundials eigen blas-devel pkg-config
   conda activate pele
   pip install git+https://github.com/pele-python/pele

If you already have GCC, G++, and GFortran, you can omit ``compilers`` from
the environment command. SUNDIALS, Eigen, and LAPACK are taken from the active
conda environment, or from the system when no conda environment is active.

For editable installs and build options, see the
`development instructions in the README <https://github.com/pele-python/pele#development>`_.

Optional features
-----------------

The GUI uses PyQt5 and PyOpenGL. With a wheel installation, add them with:

.. code-block:: console

   pip install 'pele[gui]'

The OpenGL viewer also needs the system GLU and GLUT libraries. On Ubuntu:

.. code-block:: console

   sudo apt install libglu1-mesa libglut3.12

LAMMPS potentials require the LAMMPS Python package:

.. code-block:: console

   pip install 'pele[lammps]'

Wheels use the Python LAMMPS interface. To build the compiled interface or
generate GUI forms from source, follow the
`source installation instructions in the README <https://github.com/pele-python/pele#installation>`_.

``scikit-sparse`` optionally provides sparse Cholesky for rate calculations,
and ``pymol-open-source`` provides structure viewing.
