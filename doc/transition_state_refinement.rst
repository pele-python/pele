.. _ts_refinement_description:

Transition State Refinement
===========================
This hybrid-eigenvector-following routine (a.k.a. the dimer method) finds the nearest transition state to a starting structure.
Ideally the starting point is near a transition state, but this algorithm (or a
very similar one) can be used for single ended transition state searches as
well.

A transition state is defined as a point on the energy landscape with zero
gradient and one negative eigenvalue (Higher order saddles are also possible,
but not considered here).

In order to find this point we iteratively maximize in the direction parallel
to the lowest eigenvector (the eigenvector corresponding to the lowest
eigenvalue) which minimizing in all other directions.

The algorithm consists of:

  1. :ref:`find the lowest eigenvalue <find_lowest_eigenvector_description>` and corresponding eigenvector

  2. step uphill parallel to the lowest eigenvector

  3. minimize in the space perpendicular to the lowest eigenvector

Step 3. in the above algorithm uses a standard minimization algorithm.  Step 2.
takes an eigenvector following (Newton-Raphsson-like) step uphill.  Step 1. can
be done in two ways.  If the Hessian is know, it can be diagonalized. However, even if the
Hessian is known this can be very slow.  A :ref:`better method <find_lowest_eigenvector_description>` uses an iterative scheme.

See :ref:`transition state search <transition_states_module>` for more complete documentation


