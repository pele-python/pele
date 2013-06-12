"""
.. currentmodule:: pele.angleaxis
Angle Axis Systems (`pele.angleaxis`)
============================================

This module implements routines for treating angle axis systems. An angle axis
system can be a system with aspherical potentials or a rigidid body system,
where rigid objects are treated as a collection of atoms.


Angle axis topology
++++++++++++++++++++++++++++++++++++
The basis of each Angle Axis system is a topology

.. autosummary::
   :toctree: generated/

    AATopology

"""

# dirty workaround. To not break the scripts import CoordsAdapter without moving it
from pele.utils.rbtools import CoordsAdapter
from aatopology import *
from rigidbody import *
from aamindist import *
from aasystem import *
