pele : Python Energy Landscape Explorer
+++++++++++++++++++++++++++++++++++++++

.. image:: https://github.com/pele-python/pele/actions/workflows/test.yml/badge.svg?branch=main
   :target: https://github.com/pele-python/pele/actions/workflows/test.yml
   :alt: Build Status

.. image:: https://codecov.io/gh/pele-python/pele/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/pele-python/pele
   :alt: Coverage Status

Tools for global optimization, attractor finding and energy landscape exploration.

Source code: https://github.com/pele-python/pele



.. figure:: https://raw.githubusercontent.com/pele-python/pele/main/lj38_gmin_dgraph.png

  Images: The global minimum energy structure of a 38 atom Lennard-Jones cluster.  On
  the right is a disconnectivity graph showing a visualization of the energy
  landscape.  The competing low energy basins are shown in color.

pele started as a python partial-rewriting of GMIN, OPTIM, and PATHSAMPLE: fortran
programs written by David Wales of Cambridge University and collaborators
(http://www-wales.ch.cam.ac.uk/software.html). The version started here https://github.com/pele-python/pele (documentation: http://pele-python.github.io/pele/)

The Martiniani group at New York University maintains the current version in the pele-python organization.

Description
===========
pele has tools for energy minimization, global optimization, saddle point
(transition state) search, data analysis, visualization and much more.  Some of
the algorithms implemented are:

#. Basinhopping global optimization
#. Potentials (Lennard-Jones, Morse, Hertzian, etc.) 
#. LBFGS minimization (plus other minimizers)
#. Attractor identification (Mixed Descent, CVODE)
#. Single ended saddle point search:
   - Hybrid Eigenvector Following
   - Dimer method
#. Double ended saddle point search
   - Nudged Elastic Band (NEB)
   - Doubly Nudged Elastic Band (DNEB)

#. Disconnectivity Graph visualization

#. Structure alignment algorithms

#. Thermodynamics (e.g. heat capacity) via the Harmonic Superposition Approximation

#. Transition rates analysis

Installation
============
The initial PyPI release is a source distribution and needs C/C++/Fortran
compilers, SUNDIALS, Eigen, LAPACK headers, and OpenMP. We recommend creating
a conda environment to supply these dependencies. The conda-forge package
will provide compiled binaries once its recipe has been accepted.

::

  $ conda create -n pele -c conda-forge python compilers sundials eigen blas-devel
  $ conda activate pele
  $ pip install git+https://github.com/pele-python/pele

Python 3.11 or newer is required. CI tests the latest Python release (currently 3.14) on Linux and macOS.

If the machine already has gcc, g++ and gfortran (e.g. :code:`sudo apt install gcc g++ gfortran`),
leave out :code:`compilers` for a much smaller environment. On macOS use homebrew's
gcc (Apple clang has no OpenMP or Fortran)::

  $ brew install gcc openblas
  $ CC=gcc-15 CXX=g++-15 FC=gfortran-15 pip install git+https://github.com/pele-python/pele

Optional: :code:`scikit-sparse` (sparse Cholesky for rate calculations) and
:code:`pymol-open-source` (viewing structures). The GUI (:code:`pele.gui`) uses
PyQt5 and PyOpenGL. Install them before building so Meson generates the forms::

  $ pip install PyQt5 PyOpenGL
  $ pip install --no-build-isolation '.[gui]' -Csetup-args=-Dgui=enabled

The OpenGL viewer also needs the system GLU and GLUT libraries (on Ubuntu,
:code:`sudo apt install libglu1-mesa libglut3.12`). Launch an example with
:code:`python examples/gui/ljsystem.py`.

Optional: LAMMPS potentials (:code:`pele.potentials.LAMMPSPotential` and the
:code:`pele.systems.AtomicClusterLAMMPS` system). Install the LAMMPS Python package
before building and the compiled interface is built too; without a build it falls
back to a pure Python version::

  $ pip install 'lammps[mpi]'
  $ pip install --no-build-isolation '.[lammps]' -Csetup-args=-Dlammps=enabled

See :code:`examples/gui/lammps_ljsystem.py` for a GUI example.

Development
-----------

pele is built with `meson <https://mesonbuild.com>`_ through
`meson-python <https://mesonbuild.com/meson-python/>`_. From a clone, in the same environment::

  $ pip install meson-python meson ninja cython numpy
  $ pip install --no-build-isolation -e .   # editable

Build options are meson options, passed with
:code:`-Csetup-args=...`, e.g. :code:`pip install . -Csetup-args=-Dbuildtype=debug`:

- :code:`-Dbuildtype=debug` (default :code:`release`)
- :code:`-Dcvode=disabled`: no CVODE / attractor identification; some tests will fail
- :code:`-Dnative=true`: opt in to :code:`-march=native` for local builds.
  The default is :code:`false`, for binaries that run on other machines.
- :code:`-Dgui=enabled` / :code:`disabled` (default :code:`auto`): generate the PyQt5 GUI forms; needs :code:`pyuic5` and :code:`pyrcc5`
- :code:`-Dlammps=enabled` / :code:`disabled` (default :code:`auto`): build the compiled LAMMPS potential; needs the :code:`lammps` Python package

The editable install keeps its build in :code:`build/`; pass :code:`-Cbuild-dir=...` to choose
another directory. Editable means code edits will lead to fresh rebuild for C++ code the next time 
you `import pele`. python edits will automatically reflect.

SUNDIALS, Eigen and LAPACK are taken from the active conda environment (or, without one,
from the system, e.g. :code:`sudo apt install libsundials-dev libeigen3-dev liblapack-dev liblapacke-dev`).
pele links LAPACK itself; it needs the LAPACKE package only for the header :code:`lapack.h`
(conda's :code:`blas-devel` provides both).
SUNDIALS must be built in double precision.

A :code:`CPATH`/:code:`PYTHONPATH` pointing at a pele clone takes precedence over the
installed package. If a build fails, remove the build directory before trying again::

  $ rm -rf build

Tests
=====

The project uses GitHub Actions for continuous integration (CI) testing on both Linux and macOS.
The badges at the top of this README show the current build status and code coverage.

The C++ tests use GoogleTest (the :code:`cpp_tests/gtest` submodule, or a system GoogleTest)
and need no Python packages::

  $ git submodule update --init --recursive
  $ meson setup build-tests -Dpython=disabled -Dtests=true -Dbuildtype=debug
  $ meson test -C build-tests

Add :code:`-Db_sanitize=address` for AddressSanitizer or :code:`-Db_coverage=true` for
coverage. On macOS, prefix :code:`meson setup` with :code:`CC=gcc-15 CXX=g++-15`. The
benchmarks in :code:`cpp_tests/source/benchmarks` build on request, e.g.
:code:`ninja -C build-tests cpp_tests/bench_lj`.

To run the Python tests on an installed pele::

  $ pip install pytest
  $ OMP_NUM_THREADS=1 pytest --pyargs pele

or :code:`pytest pele/` from a clone with an editable install. For coverage reporting (as in CI),
add :code:`--cov=pele --cov-report=term-missing`.

Contributing
============

See `CONTRIBUTING.md <CONTRIBUTING.md>`_.
