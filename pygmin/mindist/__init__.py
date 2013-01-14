"""
.. currentmodule:: pygmin.mindist

Tools for finding the best alignment between two structures

translational alignment
-----------------------
.. autosummary::
   :toctree: generated/

    alignCoM
    CoMToOrigin

rotational alignment
------------------------------------
The center of mass of the two structures must be at the origin

.. autosummary::
   :toctree: generated/

    alignRotation
    getAlignRotation

permutational aligment
----------------------
Finding the best alignment can be mapped onto the Assignment Problem and solved
very quickly using the Hungarian algorithm (aka munkres)

.. autosummary::
   :toctree: generated/

    findBestPermutation
    findBestPermutationRBMol

rotational + permutational alignment
------------------------------------
This cannot be solved exactly in a reasonable amount of time (The time scales as
natoms factorial).  Instead it's solved iteratively by iterating random rotation -> 
permutational alignment -> rotational alignment.  This generally produces a good alignment, but not
necessarily the optimal

.. autosummary::
   :toctree: generated/

    minPermDistStochastic
    ExactMatchCluster

A wrapper
---------
This wrapper provides the simplified interface for mindist that some 
applications require

.. autosummary::
   :toctree: generated/

    MinDistWrapper


"""
from backward_compatibility import *
from permutational_alignment import *
from exact_match import *
from minpermdist_stochastic import *
from rmsfit import *
from _minpermdist_policies import *

from _wrapper_atomiccluster import *
