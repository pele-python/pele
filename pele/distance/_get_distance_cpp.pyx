"""
# distutils: language = C++
"""
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport array_wrap_np
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

# use external c++ class
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

cdef cdef_get_distance(r1, r2, ndim, BC, dist_kwargs):
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
    BC : string
        Boundary conditions. Either 'cartesian', 'periodic' or 'lees-edwards'
    dist_kwargs : dict
        Extra arguments for the distance measurement, e.g. box (box size) or shear
    """

    assert ndim == 2 or ndim == 3, "Dimension outside the required range."
    assert BC == 'cartesian' or BC == 'periodic' or BC == 'lees-edwards', \
           "Boundary conditions undefined. Use 'cartesian', 'periodic' or 'lees-edwards'."

    cdef cppPeriodicDistance[INT2] *dist_per_2d
    cdef cppPeriodicDistance[INT3] *dist_per_3d
    cdef cppLeesEdwardsDistance[INT2] *dist_leesedwards_2d
    cdef cppLeesEdwardsDistance[INT3] *dist_leesedwards_3d
    cdef cppCartesianDistance[INT2] *dist_cart_2d
    cdef cppCartesianDistance[INT3] *dist_cart_3d

    cdef double *c_r1 = <double *>malloc(ndim * sizeof(double))
    cdef double *c_r2 = <double *>malloc(ndim * sizeof(double))
    cdef double *c_r_ij = <double *>malloc(ndim * sizeof(double))
    for i in range(ndim):
        c_r1[i] = r1[i]
        c_r2[i] = r2[i]
        c_r_ij[i] = 0

    if BC == 'periodic' or BC == 'lees-edwards':
        assert 'box' in dist_kwargs, "Required argument 'box' not defined in dist_kwargs."
        box = np.array(dist_kwargs['box']) #, dtype=float ???
        box_ = array_wrap_np(box)
        if BC == 'periodic':
            if ndim == 2:
                dist_per_2d = new cppPeriodicDistance[INT2](box_)
                dist_per_2d.get_rij(c_r_ij, c_r1, c_r2)
            else:
                dist_per_3d = new cppPeriodicDistance[INT3](box_)
                dist_per_3d.get_rij(c_r_ij, c_r1, c_r2)
        else:
            assert 'shear' in dist_kwargs, "Required argument 'shear' not defined in dist_kwargs."
            if ndim == 2:
                dist_leesedwards_2d = new cppLeesEdwardsDistance[INT2](box_, dist_kwargs['shear'])
                dist_leesedwards_2d.get_rij(c_r_ij, c_r1, c_r2)
            else:
                dist_leesedwards_3d = new cppLeesEdwardsDistance[INT3](box_, dist_kwargs['shear'])
                dist_leesedwards_3d.get_rij(c_r_ij, c_r1, c_r2)
    else:
        if ndim == 2:
            dist_cart_2d = new cppCartesianDistance[INT2]()
            dist_cart_2d.get_rij(c_r_ij, c_r1, c_r2)
        else:
            dist_cart_3d = new cppCartesianDistance[INT3]()
            dist_cart_3d.get_rij(c_r_ij, c_r1, c_r2)

    r_ij = []
    for i in range(ndim):
        r_ij.append(c_r_ij[i])
    free(c_r1)
    free(c_r2)
    free(c_r_ij)
    return r_ij

def get_distance(r1, r2, ndim, BC, dist_kwargs):
    """Python interface for the C++ get_distance implementation."""
    return cdef_get_distance(r1, r2, ndim, BC, dist_kwargs)
