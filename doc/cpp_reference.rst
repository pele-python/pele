C++ reference
=============

The C++ API reference is generated from pele's source files with Doxygen.
Breathe renders the interfaces below in Sphinx so they are searchable alongside
the Python API. The complete Doxygen reference includes public and protected
developer interfaces and functions. Private class members and file-local
helpers are excluded:

* `Class index <doxygen/annotated.html>`_
* `pele namespace and functions <doxygen/namespacepele.html>`_
* `All file-level functions <doxygen/globals_func.html>`_
* `Source files <doxygen/files.html>`_

Run ``make -C doc html`` to build the Python and C++ documentation together.
The combined site is written to ``doc/_build/html/``, with the C++ reference
under ``doxygen/``. Its navigation links back to the Python documentation.

Potential interface
-------------------

.. doxygenclass:: pele::BasePotential
   :members:
   :undoc-members:

LBFGS minimizer
---------------

.. doxygenclass:: pele::LBFGS
   :members:
   :undoc-members:

FIRE minimizer
--------------

The native modified FIRE implementation is exposed in Python as
:class:`pele.optimize.ModifiedFireCPP` and
:func:`pele.optimize.modifiedfire_cpp`. The Python implementation is
:class:`pele.optimize.Fire`, with the convenience function
:func:`pele.optimize.fire`.

.. doxygenclass:: pele::MODIFIED_FIRE
   :members:
   :undoc-members:

CVODE minimizer
---------------

.. doxygenclass:: pele::CVODEBDFOptimizer
   :members:
   :undoc-members:

Mixed descent
-------------

.. doxygenclass:: pele::GenericMixedDescent
   :members:
   :undoc-members:

Other optimizer interfaces
--------------------------

The complete Doxygen reference also documents these native optimizers:

* :cpp:class:`pele::GradientDescent`
* :cpp:class:`pele::CosineGradientDescent`
* :cpp:class:`pele::ExtendedMixedOptimizer`
* :cpp:class:`pele::MixedDescentEndOnly`

Potential interfaces
--------------------

Python constructors and methods are listed in :ref:`potentials_description`.

Pairwise potential interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: pele::PairwisePotentialInterface
   :members:
   :undoc-members:

Pairwise implementation interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Developers implementing pairwise potentials can use the direct pair loop or
the cell-list implementation. These public interfaces support the concrete
potential classes.

.. doxygenclass:: pele::SimplePairwisePotential
   :members:
   :undoc-members:
   :protected-members:

.. doxygenclass:: pele::CellListPotential
   :members:
   :undoc-members:

Potential implementations
~~~~~~~~~~~~~~~~~~~~~~~~~

The corresponding native interfaces, including their periodic and cell-list
variants, are available in the complete Doxygen reference:

* :cpp:class:`pele::LJ` and :cpp:class:`pele::LJCut`
* :cpp:class:`pele::HS_WCA` and :cpp:class:`pele::WCA`
* :cpp:class:`pele::InversePower` and :cpp:class:`pele::InversePowerHS`
* :cpp:class:`pele::InversePowerStillinger`,
  :cpp:class:`pele::InversePowerStillingerCut`, and
  :cpp:class:`pele::InversePowerStillingerCutQuad`
* :cpp:class:`pele::Frenkel` and :cpp:class:`pele::Morse`
* :cpp:class:`pele::Harmonic` and :cpp:class:`pele::RadialGaussian`
