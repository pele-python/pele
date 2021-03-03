"""
.. currentmodule:: pele.landscape
Landscape Exploration (`pele.landscape`)
============================================

This module implements routines for exploring the energy landscape.  This
primarily consists of using DoubleEndedConnect to find connected pathways
of minimum -> transition state -> minimum between two minima.

Double ended transition state search
++++++++++++++++++++++++++++++++++++
Double ended transition state searches are the main technique we use
for exploring the energy landscape.  We attempt to try to find a connected
series of minima and transition states between two end point minima.  

.. autosummary::
   :toctree: generated/

    DoubleEndedConnect

Connect manager
++++++++++++++++++++++++++++
The connect manager is a tool to organize which minima from a database 
are selected for double ended connect jobs.  

.. autosummary::
   :toctree: generated/

    ConnectManager


Other utilities
++++++++++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

    TSGraph
    database2graph
    smoothPath

Core Routines
+++++++++++++
These are some core routines used by this module.  The user probably won't need to call them,
but will want to know about them.  Parameters for these routines can be changed by passing
dictionaries to :class:`.DoubleEndedConnect`.

.. autosummary::
    :toctree: generated/
    
    LocalConnect

More core routines can be found in the documentation for the 
:ref:`transition_states <transition_states_module>` module
"""
from __future__ import absolute_import



from ._graph import *
from .local_connect import *
from .connect_min import *
#from singleended import *
from ._smooth_path import *
from .connect_manager import *

