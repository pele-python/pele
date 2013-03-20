"""
.. currentmodule:: pygmin.landscape
Landscape Exploration (`pygmin.landscape`)
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
    DoubleEndedConnectPar

Other utilities
++++++++++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

    Graph
    smoothPath

Core Routines
+++++++++++++
These are some core routines used by this module.  The user probably won't need to call them,
but will want to know about them.  Parameters for these routines can be changed by passing
dictionaries to DoubleEndedConnect

.. autosummary::
    :toctree: generated/
    
    LocalConnect
    LocalConnectPar

"""



from _graph import *
from local_connect import *
from connect_min import *
from connect_min_parallel import *
from singleended import *
from _smooth_path import *