"""
.. currentmodule:: pele.angleaxis
Angle Axis Systems (`pele.angleaxis`)
============================================

This module implements routines for treating angle axis systems. An angle axis
system can be a system with non-spherical potentials or a rigid body system,
where rigid objects are treated as a collection of atoms.

Single Rigid Body
+++++++++++++++++
The degrees of freedom of a rigid body is fully described the 3 center of mass coordinates
and the 3 elements of an angle axis rotation.  But to do calculation you need to be able
to convert from center of mass + angle axis to atomistic coordinates and back.  Furthermore
you need to know about the masses of the atoms, the tensor of gyration and more.  All of this
is stored in the RigidFragment class

.. autosummary::
   :toctree: generated/

    RigidFragment


Topology class
++++++++++++++++++++++++++++++++++++
The topology classes keep track of which types of rigid bodies, and how many of each type.
It also manages how to convert from angle axis coordinates to atomistic coordinates.

.. autosummary::
   :toctree: generated/

    RBTopology

Measurements and Transformations
++++++++++++++++++++++++++++++++
The following classes define how to perform measurements (e.g. center of mass, 
distance between structures, etc.).  They also define how to perform transformations on
a structure, e.g. how to translate a structure, how to rotate a cluster, etc.

.. autosummary::
   :toctree: generated/

    MeasureRigidBodyCluster
    TransformAngleAxisCluster

Structure alignment
++++++++++++++++++++
The following routines perform structure alignment on rigid body clusters

.. autosummary::
   :toctree: generated/

    ExactMatchAACluster
    MinPermDistAACluster


"""

# dirty workaround. To not break the scripts import CoordsAdapter without moving it
from pele.utils.rbtools import CoordsAdapter
from aatopology import AASiteType, AATopology, interpolate_angleaxis, TakestepAA
from rigidbody import RigidFragment, RBTopology, RBTopologyBulk, RigidFragmentBulk
from aamindist import TransformAngleAxisCluster, MeasureAngleAxisCluster, \
    MeasureRigidBodyCluster, ExactMatchAACluster, MinPermDistAACluster
from aasystem import AASystem, RBSystem
from _cpp_aa import RBPotentialWrapper

