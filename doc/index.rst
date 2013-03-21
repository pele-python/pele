.. pygmin documentation master file, created by
   sphinx-quickstart on Wed Aug  1 03:04:59 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pygmin's documentation!
==================================

The source code repository can be found at `<https://github.com/js850/PyGMIN>`_

The documentation is hosted at `<http://js850.github.com/PyGMIN/>`_

pygmin is a package of tools for exploring energy landscapes.  The core
routines are broken into two parts: :ref:`Basinhopping <global_optimization>`,
for finding the global minimum of an energy landscape, and for building up
databases of minima.  And :ref:`DoubleEndedConnect <landscape_module>`, for
finding minimum energy paths on the energy landscape between two minima.  This
means paths that go through the geometric transition states.

Note that we use the language "energy landscape" because it's the language most
natural for our fields of physics and chemistry, but most of these tools are
equally applicable for working with any smooth scalar function in N dimensions.

Tutorials
-----------
.. toctree::
   :maxdepth: 3

   tutorial_potential
   system_class_tutorial


Reference
---------

.. toctree::
   :maxdepth: 2

   global_optimization

.. toctree::
   :maxdepth: 3

   landscape

.. toctree::
   :maxdepth: 2

   disconnectivity_graph

Modules
+++++++
.. toctree::
   :maxdepth: 1

   system_class
   quenching
   potentials
   database
   structure_alignment
   step_taking
   transition_states_module
   landscape_module
   monte_carlo
   accept_tests
   utils
   gui

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

