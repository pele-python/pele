"""
.. currentmodule:: pygmin.mindist

Tools for finding the best alignment between two structures

rotational alignment
------------------------------------
Find the optimal rotation for two sets of coordinates (i.e. perform a rms fit)

.. autosummary::
   :toctree: generated/

    findrotation
    findrotation_kabsch
    findrotation_kearsley

permutational aligment
----------------------
Finding the best alignment can be mapped onto the Assignment Problem and solved
very quickly using the Hungarian algorithm (aka munkres) or shortest augmenting path

.. autosummary::
   :toctree: generated/

    optimize_permutations
    find_best_permutation
    find_permutations_munkres
    find_permutations_hungarian
    find_permutations_OPTIM
    

rotational + permutational alignment
------------------------------------
This cannot be solved exactly in a reasonable amount of time (The time scales as
natoms factorial).  Instead it's solved iteratively by iterating random rotation -> 
permutational alignment -> rotational alignment.  This generally produces a good alignment, but not
necessarily the optimal

.. autosummary::
   :toctree: generated/

    MinPermDistCluster
    ExactMatchCluster
    StandardClusterAlignment

For atomic cluster, specialized wrapper exist.

.. autosummary::
   :toctree: generated/
    
    MinPermDistAtomicCluster
    ExactMatchAtomicCluster

See the angleaxis module for angleaxis minpermdist routines

Periodic Boundary Conditions
----------------------------
This is generally a much harder problem than those discussed above.  Currently
we have no general mindist routine, but we do have a test to check if the
strructures are identical

.. autosummary::
   :toctree: generated/
    
    ExactMatchPeriodic
    
Customizing minpermdist - minpermdist policies
----------------------------------------------
.. autosummary::
   :toctree: generated/

    TransformPolicy
    MeasurePolicy
    
Utilities
---------
.. autosummary::
   :toctree: generated/

    PointGroupOrderCluster


OBSOLETE: translational alignment
-----------------------
.. autosummary::
   :toctree: generated/

    alignCoM
    CoMToOrigin

"""
from backward_compatibility import *
from permutational_alignment import *
from exact_match import *
from minpermdist_stochastic import *
from rmsfit import *
from _minpermdist_policies import *
from periodic_exact_match import ExactMatchPeriodic
from _pointgrouporder import *
from _wrapper_atomiccluster import *
