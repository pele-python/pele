.. _landscape_description:

Energy Landscape Exploration
============================

This is an attempt at describing the methods used in the double ended connect
algorithm.  For instructions on how to use it see the :ref:`module
<landscape_module>` documention and the documentation of functions therin.
Also, there are examples in the examples/ folder describing how to use these
routines.


Double ended transition state search
++++++++++++++++++++++++++++++++++++

.. toctree::
  :maxdepth: 1

  local_connect
  structure_alignment
  
Double ended transition state searches are the main technique we use
for exploring the energy landscape.  We attempt to try to find a connected
series of transition states between two minima.  The end result might look like::

                             TS2
                TS1        /    \
              /    \      /      \
             /       Min1         \
    MinStart                       \
                                    MinEnd

The global algorithm looks like this::

    while MinStart and MinEnd are not connected:
        min1, min2 = get_next_minima_pair()
        local_connect(min1, min2)

The algorithm has two parts: get_next_minima_pair, and local_connect (psuedo-code names).
get_next_minima_pair selects a pair of known, but, unconnected minima (min1 and min2) to try to connect.
The algorithm which chooses this minima pair is designed to find the optimial path
as quickly as possible between MinStart and
MinEnd.
:ref:`Local connect <local_connect_description>` 
does a nudged elastic band run between min1 and min2 in order to get transition
state candidates.  It then refines the transition state candidates to find actual transition
states to the required tolerance.  It then falls off each side of the transition state
and minimized in order to find the pair of minima which the transition state connects.



See :ref:`landscape <landscape_module>` for more complete documentation
