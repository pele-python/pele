.. _find_lowest_eigenvector_description:

Finding the Lowest Eigenvector
------------------------------
This describes an iterative scheme to find the lowest eigenvalue and
corresponding eigenvector.  The scheme is simply an optimization problem:  find
the vector along which the second derivative is minimal.  The second derivative
is approximated using three points along the vector.  

Documentation
--------------
.. autofunction:: pygmin.transition_states.findLowestEigenVector

