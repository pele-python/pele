.. _local_connect_description:
Local Connect Algorithm
=======================

.. toctree::
  :maxdepth: 1

  NEB
  transition_state_refinement
  transition_state_refinenement

Local connect has three steps.  

1. Do a :ref:`nudged elastic band (NEB) <neb_description>` run between min1 and min2 in order to get
   transition state candidates.  

#. :ref:`Refine the transition state <ts_refinement_description>` candidates to find actual transition states to
   the required tolerance.  

#. Falls off both sides of the transition state and minimize in order to find
   the pair of minima which the transition state connects.

See :ref:`landscape <landscape_module>` for more complete documentation


