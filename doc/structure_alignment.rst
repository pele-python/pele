.. _structure_alignment

Structure Alignment
===================
(a.k.a. minimum distance routines, mindist, minpermdist, etc.) When trying to
find path between two minima it is importatnt to ensure that the starting
points are as close together as possible given the symmetries of the system.
Common symmetries that need to be accounted for are global rotational
invariance, global inversion symmetry, global translational invariance, and
permutational invariance.

Translation and inversion symmetries are trivial to deal with.
Given *either* rotational or permutational symmetries it is trivial to find the
optimum alignment.  When you have both, a solution cannot be found analytically
without going through all N-factorial options, a ludicrously slow option.
Instead we solve it approximately and stochastically.


Documentation
-------------
.. automodule:: pygmin.mindist
