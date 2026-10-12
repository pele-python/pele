"""
.. currentmodule:: pele.optimize

Optimizers (`pele.optimize`)
================================
This module contains all of the localoptimizers available in pele.  There are so
many available for testing purposes and because sometimes different optimizers are
more appropriate in different circumstances.  These optimizers are all local.  That
is they stop at a local minimum.  For global optimization see :ref:`global
optimization <global_optimization>`.

In our experience, we have found that the best minimizers are lbfgs and fire.
lbfgs seems to be by far the fastest, however there are some circumstances
where the better stability of fire might be useful. CVODE exists if you want guarantees
that the endpoint is the one connected by the steepest descent equation (see https://arxiv.org/abs/2409.12113).
Mixed Descent is also designed to increase the probability that the path is connected by steepest descent,
(but is currently WIP because it only improves performance over CVODE for LJ clusters)
LBFGS can also fail in non-Hamiltonian situations where you have
forces, but no potential function.  This is the case in the Nudged elastic band,
where the NEB force doesn't correspond to any NEB energy.

.. note::
    we use the words "optimize", "quench", and "minimize" interchangeably


All of the optimizers that are functions in this package have the same form::

    res = optimizer(coords, potential, **kwargs)

where ``coords`` is the starting structure, ``potential`` supplies the energy and
gradient, and ``kwargs`` contains optional parameters.
The return value is a :class:`Result` object,
which is simply a dictionary where `__getattr__` is a wrapper for `__getitem__`.  So the
coords can be accessed as res.coords or as res["coords"].

.. autosummary::
   :toctree: generated/

   Result

We have tried to make the minimizers as consistent as possible, but it is not always possible.
Some of the common parameters most of them accept are:

1. nsteps : number of iterations
#. tol : tolerance criterion for the rms gradient
#. iprint : how often to print status information

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


C++ LBFGS
~~~~~~~~~
:func:`lbfgs_cpp` constructs :class:`LBFGS_CPP` and returns its result.
The system class uses this implementation for its default local minimizer.

.. autosummary::
   :toctree: generated/

   LBFGS_CPP
   lbfgs_cpp

Python and Fortran LBFGS
~~~~~~~~~~~~~~~~~~~~~~~~
:class:`LBFGS` is the Python implementation; :class:`MYLBFGS` uses a Fortran
step update. Their function wrappers return the optimization result.

.. autosummary::
   :toctree: generated/

   LBFGS
   lbfgs_py
   MYLBFGS
   mylbfgs

FIRE
----

Python FIRE
~~~~~~~~~~~
:func:`fire` runs the Python :class:`Fire` implementation.

.. autosummary::
    :toctree: generated/

    Fire
    fire

Native modified FIRE
~~~~~~~~~~~~~~~~~~~~
:func:`modifiedfire_cpp` runs :class:`ModifiedFireCPP`, the Python interface
to the C++ ``pele::MODIFIED_FIRE`` implementation::

    from pele.optimize import modifiedfire_cpp
    result = modifiedfire_cpp(coords, potential)

To retain the optimizer object between calls, construct
``ModifiedFireCPP(coords, potential, **kwargs)`` and call its ``run()`` method.

.. autosummary::
    :toctree: generated/

    ModifiedFireCPP
    modifiedfire_cpp

Native gradient descent
-----------------------
These C++ optimizer classes take the potential before the initial coordinates.
Call ``run()`` on the optimizer to obtain its result.

.. autosummary::
    :toctree: generated/

    GradientDescent_CPP
    CosineGradientDescent

CVODE and mixed descent
-----------------------
These native optimizers require a build with SUNDIALS and the ``cvode`` feature
enabled. Their constructors take the potential before the initial coordinates.
:class:`GenericMixedDescent` also takes the non-convex and convex optimizer
objects, a tolerance, and a translation offset.

.. autosummary::
    :toctree: generated/

    CVODEBDFOptimizer
    MXDOptimizer
    ExtendedMixedOptimizer
    GenericMixedDescent

Native optimizer options
------------------------
:class:`StopCriterionType` selects the stopping criterion for native LBFGS,
modified FIRE, and extended mixed descent. :class:`HessianType` selects the
Hessian implementation used by :class:`CVODEBDFOptimizer`.

.. autosummary::
    :toctree: generated/

    StopCriterionType
    HessianType

SciPy optimizers
----------------
These functions wrap the corresponding optimizers in SciPy.

.. autosummary::
    :toctree: generated/

    lbfgs_scipy
    bfgs_scipy
    cg

Python steepest descent
-----------------------
This implementation is intended for testing and is slow.

.. autosummary::
    :toctree: generated/

    steepest_descent


"""

from .result import *
from ._lbfgs_py import *
from ._mylbfgs import *
from ._fire import *
from ._gradient_descent_cpp import GradientDescent_CPP
from ._modified_fire_cpp import ModifiedFireCPP
from ._lbfgs_cpp import LBFGS_CPP, StopCriterionType
from ._cosine_gradient_descent_cpp import (
    CosineGradientDescent_CPP as CosineGradientDescent,
)

try:
    from .cvode_opt import (
        CVODEBDFOptimizer_CPP as CVODEBDFOptimizer,
        HessianType,
    )
    from ._mxd_end_only import MXDOptimizer_CPP as MXDOptimizer
    from .extended_mixed_descent import (
        ExtendedMixedOptimizer_CPP as ExtendedMixedOptimizer,
    )
    from .generic_mixed_descent import (
        GenericMixedDescent_CPP as GenericMixedDescent,
    )
except ImportError:
    print("Warning: CVODE based optimizers not available")
    pass

from ._quench import *
