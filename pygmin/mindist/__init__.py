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

.. autosummary::
   :toctree: generated/

    findBestPermutation
    findBestPermutationRBMol

rotational + permutational alignment
------------------------------------
This cannot be solved exactly in a reasonable amount of time (The time scales as
natoms factorial).  Instead it's treated as a global optimization problem and solved
stochastically using basinhoping.  This generally produces a good alignment, but not
necessarily the optimal

.. autosummary::
   :toctree: generated/

    minPermDistStochastic


"""

from minpermdist import *
from minpermdist_stochastic import *
from minpermdist_rbmol import *
from mindistutils import *
from rmsfit import *