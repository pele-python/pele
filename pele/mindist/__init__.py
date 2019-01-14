"""
.. currentmodule:: pele.mindist

Structure Alignment (`pele.mindist`)
=========================================
Tools for finding the best alignment between two structures.
(a.k.a. minimum distance routines, mindist, minpermdist, etc.) When trying to
find path between two minima it is important to ensure that the starting
points are as close together as possible given the symmetries of the system.
Common symmetries that need to be accounted for are global translational 
invariance, global rotational invariance, global inversion symmetry, and
permutational invariance.

Translation and inversion symmetries are trivial to deal with.
Given *either* rotational or permutational symmetries it is trivial to find the
optimum alignment.  When you have both, a solution cannot be found analytically
without going through all N-factorial options, a ludicrously slow option.
Instead we solve it approximately and stochastically.

All `mindist` routines that are used throughout pele should have the same signature::

    distance x1new, x2new = compute_mindist(x1, x2)

Where it is important that `x1new` is the *same* as `x1`.  `x2new` is `x2` after being
put into best alignment with `x1`.  This signature is expected wherever `mindist` routines are used
and is expected in :meth:`.BaseSystem.get_mindist`.
Often `compute_mindist` is a class that has overloaded `__call__` as is the case with, e.g.::

    compute_mindist = MinPermDistAtomicCluster(permlist=[range(natoms)])
    distance x1new, x2new = compute_mindist(x1, x2)

    



Rotational alignment
------------------------------------
Find the optimal rotation for two sets of coordinates (i.e. perform a rms fit)

.. autosummary::
   :toctree: generated/

    findrotation
    findrotation_kabsch
    findrotation_kearsley

Permutational aligment
----------------------
Finding the best alignment can be mapped onto the Assignment Problem and solved
very quickly using the Hungarian (shortest augmenting path) algorithm.

.. autosummary::
   :toctree: generated/

    optimize_permutations
    find_best_permutation
    

Rotational + permutational alignment
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
structures are identical

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

"""
from __future__ import absolute_import
from .backward_compatibility import CoMToOrigin
from .permutational_alignment import *
from .exact_match import StandardClusterAlignment, ExactMatchCluster
from .minpermdist_stochastic import MinPermDistCluster
from .rmsfit import *
from ._minpermdist_policies import *
from .periodic_exact_match import ExactMatchPeriodic, MeasurePeriodic, TransformPeriodic
from ._pointgrouporder import PointGroupOrderCluster
from ._wrapper_atomiccluster import MinPermDistAtomicCluster, ExactMatchAtomicCluster

