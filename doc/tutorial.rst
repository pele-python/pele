pygmin Turorial
===============

pygmin is a package for global optimization and energy landscape analysis.
For global optimization the preferred method is to use basinhopping.
basinhopping a stochastic search algorithm whos structure is very similar to
the Monte Carlo algorithm but does a local minimization (or quench) after every
step and uses the quenched energy in the accept/reject criterion.

By energy landscape analysis we generally mean finding geometrical transition
states between known minima.  If we build up a network of minima and transition
states we can have a very good picture of the thermodynamic and even dynamic
properties of the system.  A general work cycle often looks like

1. use :ref:`Basinhopping <global_optimization>` to find the global minimum and other important minima

2. use :ref:`DoubleEndedConnect <landscape_description>` to find the connected pathways of minima and
   transition states between these known minima

3. post anysis might include plotting the :ref:`disconnectivity graph <disconnectivity_graph>`



.. .. include:: tutorial_potential.rst

.. .. include:: system_class_tutorial.rst
