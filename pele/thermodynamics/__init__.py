"""
.. currentmodule:: pele.thermodynamics
Thermodynamics (`pele.thermodynamics`)
=========================================

This module implements routines to calculate thermodynamic properties within
the harmonic apprixiamtion and superposition principle. 

Normal mode analysis
++++++++++++++++++++
Most of the calculations (free energy, heat capacity curves, ...) are based on
normal mode analysis

.. autosummary::
   :toctree: generated/

    normalmode_frequencies
    logproduct_freq2

Heat Capacity
+++++++++++++
Heat capacity and other quantities can be calculated from 
a database of minima once the normal mode frequencies and the
point group order are known.  

.. autosummary::
   :toctree: generated/

    dos_to_cv
    minima_to_cv

Utilities
---------
These are functions which you may find useful.

.. autosummary::
   :toctree: generated/

    get_thermodynamic_information

    

"""
from __future__ import absolute_import

from ._normalmodes import *
from .heat_capacity import *
from ._utils import *

