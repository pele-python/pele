"""
.. currentmodule:: pygmin.systems

System Class
============

The system class is a convenience wrapper for easily defining all the necessary
information about a given system in one place.  It also makes using the tools in
pygmin very simple.


BaseSystem Class
----------------

the BaseSystem defines the base class from which all other system
classes are derived.  It lists the optional and required functions for defining
a system class.  It also defines several functions (e.g. get_basinhopping, get_double_ended_connect)
which uses other functions to create high level objects.  Depending on what calculations
you want to perform, different functions are required

basinhopping::

        get_potential : required
        get_takestep : optional
        get_random_configuration : optional
        get_compare_exact : optional

landscape exploration and transition state searches::

        get_potential : required
        get_mindist : required
        get_orthogonalize_to_zero_eigenvectors : required
        get_compare_exact : optional, recommended
        get_random_configuration : optional
        
See :ref:`Potentials <potentials_description>` for more information about how to implement get_potential().
See :ref:`Structure Alignment <structure_alignment>` for more information about how to implement get_mindist().
See :ref:`Global Optimization <global_optimization>` for more information about how to implement get_takestep().

 

.. autosummary::
    :toctree: generated
    
    BaseSystem


classes of systems
------------------
These are groups of system types that have common features (e.g. the same 
symmetries).  Thus some of the required functions can be defined here 
and inherited by specific systems.

.. autosummary::
    :toctree: generated
    
    AtomicCluster

system list
-----------
Here we list the existing systems.

.. autosummary::
    :toctree: generated
    
    LJCluster
    BLJCluster

Parameter Class
---------------
The parameter class is what the system class uses for
holding and maintaining global defaults for all the many adjustable
parameters in pygmin.  

If, in your new system class you want to change the default number of
double ended connect iterations, you can do it by adding this line
to __init__()::

    self.params.double_ended_connect.niter = 50

then, when you create the double ended connect object it will use
this new value::
    
    connect = mysys.get_double_ended_connect(min1, min2, database)


The top level functions often use multiple levels of algorithms.  E.g. 
DoubleEndedConnect calls LocalConnect which calls FindTransitionState.  
If you want to change a parameter in one of these low level algorithm
you must use the parameter tree.
the parameter tree looks like this::

    -params
    -----double_ended_connect
    ---------local_connect_params
    -------------NEBparams
    -----------------NEBquenchParams
    -------------tsSearchParams
    -----------------lowestEigenvectorQuenchParams
    -----------------tangentSpaceQuenchParams

            
All of the these listed above are dictionaries that are passed to the appropriate functions or classes.  
For simplification, they 
can be accessed either by keyword or attribute.  The maximum uphill step in the transition state search can
be modified in either of these two equivalent ways::

    mysys.params.double_ended_connect.local_connect_params.tsSearchParams.max_uphill_step = 0.2
    mysys.params["double_ended_connect"]["local_connect_params"]["tsSearchParams"]["max_uphill_step"] = 0.2

The logic is that the system class passes the dictionary `double_ended_connect`
as kwargs (key word arguments) to `DoubleEndedConnect`.  One of those keyword
argemnts is `local_connect_params`, which is then passed by
`DoubleEndedConnect` as kwargs to `LocalConnect`.  `LocalConnect` then passes
`tsSearchParams` to `FindTransitionState` and so on.  Each parameter dictionary
holds keyword arguments for the following function

===============================  =============================================
parameter dictionary             passed as parameters to
===============================  =============================================
`double_ended_connect`           :ref:`DoubleEndedConnect <landscape_module>`
`local_connect_params`           :ref:`LocalConnect <landscape_module>`
`NEBparams`                      :ref:`create_NEB <transition_states_module>`
`NEBquenchParams`                the :ref:`optimizer <optimize_module>` called by :ref:`NEB <transition_states_module>`
`tsSearchParams`                 :ref:`FindTransitionState <transition_states_module>`
`lowestEigenvectorQuenchParams`  the :ref:`optimizer <optimize_module>` called by :ref:`findLowestEigenVector <transition_states_module>`
`tangentSpaceQuenchParams`       the :ref:`optimizer <optimize_module>` called by :ref:`FindTransitionState <transition_states_module>` 
                                   for the tangent space quench
===============================  =============================================




Note that the Parameters class doesn't hold the pygmin default values for
each algorithm.  These are defined in the algorithms themselves.  Instead,
the Parameters class keeps track only of what has been overridden.

.. autosummary::
    :toctree: generated
    
    BaseParameters
    Parameters
    
For a translation between an OPTIM odata file and the pygmin Parameter tree, see
:ref:`here <optim2params>`

"""


from basesystem import *
from cluster import *
from ljcluster import *
from bljcluster import *

