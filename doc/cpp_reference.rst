C++ reference
=============

The C++ API reference is generated from pele's source files with Doxygen.
Breathe renders the interfaces below in Sphinx so they are searchable alongside
the Python API. The complete Doxygen reference includes all extracted classes
and functions:

* `Class index <doxygen/annotated.html>`_
* `pele namespace and functions <doxygen/namespacepele.html>`_

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
