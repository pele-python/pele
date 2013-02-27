"""
.. currentmodule:: pygmin.thermodynamics
Thermodynamics
=====================

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
"""

from _normalmodes import *