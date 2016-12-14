"""
# distutils: language = C++
"""
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport array_wrap_np
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

# cython has no support for integer template argument.  This is a hack to get around it
# https://groups.google.com/forum/#!topic/cython-users/xAZxdCFw6Xs
# Basically you fool cython into thinking INT2 is the type integer,
# but in the generated c++ code you use 2 instead.
# The cython code MyClass[INT2] will create c++ code MyClass<2>.
cdef extern from *:
    ctypedef int INT2 "2"    # a fake type
    ctypedef int INT3 "3"    # a fake type

# use external c++ classes
cdef extern from "pele/distance.h" namespace "pele":
    cdef cppclass cppCartesianDistance "pele::cartesian_distance"[ndim]:
        cppCartesianDistance() except +
        void get_rij(double *, double *, double *)
    cdef cppclass cppPeriodicDistance "pele::periodic_distance"[ndim]:
        cppPeriodicDistance(_pele.Array[double] box) except +
        void get_rij(double *, double *, double *)
    cdef cppclass cppLeesEdwardsDistance "pele::leesedwards_distance"[ndim]:
        cppLeesEdwardsDistance(_pele.Array[double] box, double shear) except +
        void get_rij(double *, double *, double *)

cpdef get_distance(np.ndarray[double] r1, np.ndarray[double] r2, int ndim, str method, box=None, double shear=0.0):
    """
    Define the Python interface to the C++ distance implementation.

    Parameters
    ----------
    r1 : [float]
        Position of the first particle
    r2 : [float]
        Position of the second particle
    ndim : int
        Number of dimensions
    method : string
        Distance measurement method / boundary conditions. Either 'cartesian', 'periodic' or 'lees-edwards'
    box : np.array[float], optional
        Box size (for periodic and Lees-Edwards distance measure)
    shear : float, optional
        Amount of shear (for Lees-Edwards distance measure)
    """

    # Assert that the input parameters are right
    assert ndim == 2 or ndim == 3, "Dimension outside the required range."
    assert method == 'cartesian' or method == 'periodic' or method == 'lees-edwards', \
           "Distance measurement method undefined. Use 'cartesian', 'periodic' or 'lees-edwards'."

    # Define pointers for all distance measures
    # (otherwise this would be clumsily handled by Cython, which would call
    #  the empty constructor before constructing the objects properly.)
    cdef cppPeriodicDistance[INT2] *dist_per_2d
    cdef cppPeriodicDistance[INT3] *dist_per_3d
    cdef cppLeesEdwardsDistance[INT2] *dist_leesedwards_2d
    cdef cppLeesEdwardsDistance[INT3] *dist_leesedwards_3d
    cdef cppCartesianDistance[INT2] *dist_cart_2d
    cdef cppCartesianDistance[INT3] *dist_cart_3d

    # Define box in C
    cdef _pele.Array[double] c_box

    # Initialize data arrays in C
    cdef double *c_r1 = <double *>malloc(ndim * sizeof(double))
    cdef double *c_r2 = <double *>malloc(ndim * sizeof(double))
    cdef double *c_r_ij = <double *>malloc(ndim * sizeof(double))
    for i in range(ndim):
        c_r1[i] = r1[i]
        c_r2[i] = r2[i]

    # Calculate the distance
    if method == 'periodic' or method == 'lees-edwards':

        # Get box size from the input parameters
        assert box is not None, "Required argument 'box' not defined."
        c_box = array_wrap_np(box)

        if method == 'periodic':
            if ndim == 2:
                dist_per_2d = new cppPeriodicDistance[INT2](c_box)
                dist_per_2d.get_rij(c_r_ij, c_r1, c_r2)
            else:
                dist_per_3d = new cppPeriodicDistance[INT3](c_box)
                dist_per_3d.get_rij(c_r_ij, c_r1, c_r2)
        else:
            if ndim == 2:
                dist_leesedwards_2d = new cppLeesEdwardsDistance[INT2](c_box, shear)
                dist_leesedwards_2d.get_rij(c_r_ij, c_r1, c_r2)
            else:
                dist_leesedwards_3d = new cppLeesEdwardsDistance[INT3](c_box, shear)
                dist_leesedwards_3d.get_rij(c_r_ij, c_r1, c_r2)
    else:
        if ndim == 2:
            dist_cart_2d = new cppCartesianDistance[INT2]()
            dist_cart_2d.get_rij(c_r_ij, c_r1, c_r2)
        else:
            dist_cart_3d = new cppCartesianDistance[INT3]()
            dist_cart_3d.get_rij(c_r_ij, c_r1, c_r2)

    # Copy results into Python object
    r_ij = np.empty(ndim)
    for i in xrange(ndim):
        r_ij[i] = c_r_ij[i]

    # Free memory
    free(c_r1)
    free(c_r2)
    free(c_r_ij)

    return r_ij
