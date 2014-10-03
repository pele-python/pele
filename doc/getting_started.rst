Getting started
===============

pele is a package for global optimization and landscape analysis.
The general sequence of steps when using these methods is

1. Find the global minimum and the build up a database of other important minima

   - We use :ref:`Basinhopping <global_optimization>` to do the global optimization.

#. Find connections between the minima.  

   - This means finding the saddle point separating two minima such that the
     minimum energy path between the two minima crosses through the saddle
     point.

   - We use :ref:`DoubleEndedConnect <landscape_description>` which is a combination
     of the Nudged Elastic Band (:class:`.NEBDriver`) and Hybred Eigenvector following (:class:`.FindTransitionState`)

#. Landscape analysis.

   - Visualize the landscape with a :ref:`disconnectivity graph <disconnectivity_graph>`.
   - Compute :ref:`rates <rates_module>` between minima or groups of minima



.. automodule:: examples.getting_started.example1_minimization
