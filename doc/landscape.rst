Energy Landscape Exploration
============================


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
get_next_minima_pair selects a pair of known minima that are unconnected with the ultimate
goal of connecting MinStart and MinEnd.  See here for a description of this routine.
:ref:`Local connect <local_connect_description>` 
does a nudged elastic band run between min1 and min2 in order to get transition
state candidates.  It then refines the transition state candidates to find actual transition
states to the required tolerance.  It then falls off each side of the transition state
and minimized in order to find the pair of minima which the transition state connects.

We start by trying to connect MinStart and MinEnd, and as new transition states and minima are added
As we find new tran
The main algorithm 


.. currentmodule:: pygmin.landscape

.. autoclass:: DoubleEndedConnect
