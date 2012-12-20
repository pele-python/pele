""" 
.. currentmodule:: pygmin.optimize

Optimizers
==========
This module contains all of the optimizers available in pygmin.  There are so 
many available for testing purposes and because sometimes different optimizers are
more appropriate in different circumstances.  These optimizers are all local.  That
is they find the nearest local minimum.  For global optimization see :ref:`global
optimization <global_optimization>`.

In our experience, we have found that the best minizers are lbfgs and fire.  
lbfgs seems to be by far the fastest, however there are some circumstances
where the use of lbfgs is inapropriate, and we suggest fire instead.  These cases
include where it is important to clearly define the boundary between basins of 
attraction.  LBFGS can also fail in non-Hamiltonian situations where you have
gradient information, but no potential function.  This is the case in the Nudged 
elastic band, where the NEB force doesn't correspond to any NEB energy.

.. note::
    we use the words "optimize", "quench", and "minimize" interchangeably


All of the optimizers that are functions in this package have the same form::
    
    res = optimizer(coords, getEnergyGradient, **kwargs)

where `coords` is the starting structure for the minimaztion, getEnergyGradient is a function 
which return the energy and gradient, and kwargs are a collection of optional parameters.
We have tried to make the minimizers as consistent as possible, but it is not always possible.
Some of the common parameters most of them accept are::

1. nsteps : number of iterations
#. tol : tolerance criterion for the rms gradient
#. iprint : how often to print status informatation 

LBFGS routines
--------------
limited-memory BFGS (Broyden-Fletcher-Goldfarb-Shanno) routine

http://en.wikipedia.org/wiki/Limited-memory_BFGS

These routines (excluding the scipy version) are slightly different from
most lbfgs routines in that the don't use a line search.  Instead, both
the step size and direction returned by the lbfgs routine are accepted subject
to a constraint on the energy change.  If the energy rises more than a given amount then the 
step size is reduce until the condition is satisfied.  Note: this is what makes
lbfgs potentially fail with non-Hamiltonian systems.

.. autosummary::
   :toctree: generated/
   
   LBFGS
   lbfgs_py
   MYLBFGS
   mylbfgs
   lbfgs_scipy

Fire
----
.. autosummary::
    :toctree: generated/
    
    Fire
    fire


Other routines
---------------
most of these are simply wrappers to the
optimizers available in scipy.optimize.
These are not used very often and may be buggy.  

.. autosummary::
    :toctree: generated/
    
    bfgs
    cg
    fmin
    steepest_descent


"""

from result import *
from _lbfgs_py import *
from _mylbfgs import *
from _fire import *
from _quench import *