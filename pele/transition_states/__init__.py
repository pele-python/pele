"""
.. currentmodule:: pele.transition_states
Transition States (`pele.transition_states`)
================================================
This module contains functions and classes related to local transition state searches.
Most of these algorithms will not need to be directly called by the user.  However
it is important to know how these work because they form some of the core routines of
:ref:`landscape exploration <landscape_module>`.  
 
 
Local transition state search
+++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/
    
    FindTransitionState
    findTransitionState

Lowest eigenvalue search
++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

    findLowestEigenVector

(Doubly-) Nudged Elastic Band
+++++++++++++++++++
The following functions implement the Nudged Elastic Band and Doubly-Nudged Elastic Band method.
The user should interact via the driver class:
 
.. autosummary::
    :toctree: generated/
    
    NEBDriver

In the backend, the work is done by the following functions. 

.. autosummary::
    :toctree: generated/
    
    NEB
    InterpolatedPath
    InterpolatedPathDensity

Orthogonalize to zero eigenvectors
++++++++++++++++++++++++++++++++++
These functions make a vector orthogonal to known zero eigenvectors (eigenvectors with a zero eigenvalue).  Typically these
correspond to known symmetries of the system like translational invariance, rotational invariance.
Frozen degrees of freedom also contribute zero eigenvecotrs

.. autosummary::
    :toctree: generated/

    orthogopt
    orthogopt_translation_only
    zeroEV_translation
    zeroEV_rotation
    zeroEV_cluster
    gramm_schmidt

Finding minima on either side of a transition state
---------------------------------------------------
When a transition state is found we step off either side of the transition state to find
the minima which the transition state connects.  This routine controls that process

.. autosummary::
    :toctree: generated/

    minima_from_ts


"""
from __future__ import absolute_import

from ._zeroev import *
from ._orthogopt import *
from ._interpolate import *
from ._NEB import *
from ._find_lowest_eig import *
from ._transition_state_refinement import *
from ._tstools import *
from ._nebdriver import *
#from _generalized_dimer import GeneralizedDimer

