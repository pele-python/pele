"""
.. currentmodule:: pele.rates
Rates (`pele.rates`)
====================

This module contains tools to transition rates from one part of the graph to
another.  The rates are computed using the New Graph Transformation (NGT)
method of David Wales (J. Chem. Phys. 2009) http://dx.doi.org/10.1063/1.3133782

The method uses a graph renormalization method (renormalization in the sense of
renormalization group theory) to compute exact Kinetic Monte Carlo rates and
first passage probabilities from a reactant group A and the product group B.
The intermediate nodes (those not in the reactant or product groups) are
removed iteratively. Upon removing a node the waiting times and transition
probabilities of the neighbors are updated.

.. autosummary::
   :toctree: generated/
    
    GraphReduction

"""
from _rates import *
