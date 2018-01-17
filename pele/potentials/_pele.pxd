cimport numpy as np
import numpy as np

from ctypes import c_size_t as size_t

#===============================================================================
# shared pointer
#===============================================================================
cdef extern from "<memory>" namespace "std":
    cdef cppclass shared_ptr[T]:
        shared_ptr() except+
        shared_ptr(T*) except+
        T* get() except+
        # T& operator*() # doesn't do anything
        # Note: operator->, operator= are not supported

#===============================================================================
# pele::Array
#===============================================================================
cdef extern from "pele/array.h" namespace "pele":
    cdef cppclass Array[dtype] :
        Array() except +
        Array(size_t) except +
        Array(dtype*, size_t n) except +
        size_t size() except +
        Array[dtype] copy() except +
        dtype *data() except +
        dtype & operator[](size_t) except +

#===============================================================================
# pele::BasePotential
#===============================================================================
cdef extern from "pele/base_potential.h" namespace "pele":
    cdef cppclass  cBasePotential "pele::BasePotential":
        cBasePotential() except +
        double get_energy(Array[double] &x) except +
        double get_energy_gradient(Array[double] &x, Array[double] &grad) except +
        double get_energy_gradient_hessian(Array[double] &x, Array[double] &g, Array[double] &hess) except +
        void get_hessian(Array[double] &x, Array[double] &hess) except +
        void numerical_gradient(Array[double] &x, Array[double] &grad, double eps) except +
        void numerical_hessian(Array[double] &x, Array[double] &hess, double eps) except +

#cdef extern from "potentialfunction.h" namespace "pele":
#    cdef cppclass  cPotentialFunction "pele::PotentialFunction":
#        cPotentialFunction(
#            double (*energy)(Array[double] x, void *userdata) except *,
#            double (*energy_gradient)(Array[double] x, Array[double] grad, void *userdata) except *,
#            void *userdata) except +

#===============================================================================
# cython BasePotential
#===============================================================================
cdef class BasePotential:
    cdef shared_ptr[cBasePotential] thisptr      # hold a C++ instance which we're wrapping

#===============================================================================
# pele::CombinedPotential
#===============================================================================
cdef extern from "pele/combine_potentials.h" namespace "pele":
    cdef cppclass  cCombinedPotential "pele::CombinedPotential":
        cCombinedPotential() except +
        double get_energy(Array[double] &x) except +
        double get_energy_gradient(Array[double] &x, Array[double] &grad) except +
        void add_potential(shared_ptr[cBasePotential] potential) except +


#===============================================================================
# conversion routines between numpy arrays and pele::Array
#===============================================================================
cdef inline Array[double] array_wrap_np(np.ndarray[double] v) except *:
    """return a pele Array which wraps the data in a numpy array
    
    Notes
    -----
    we must be careful that we only wrap the existing data
    """
    if not v.flags["FORC"]:
        raise ValueError("the numpy array is not c-contiguous.  copy it into a contiguous format before wrapping with pele::Array")
    return Array[double](<double *> v.data, v.size)

cdef inline np.ndarray[double, ndim=1] pele_array_to_np(Array[double] v):
    """copy the data in a pele::Array into a new numpy array
    """
    cdef int i
    cdef int N = v.size()
    cdef np.ndarray[double, ndim=1] vnew = np.zeros(N)
    for i in xrange(N):
        vnew[i] = v[i]
    return vnew
    
cdef inline Array[size_t] array_wrap_np_size_t(np.ndarray[size_t] v) except *:
    """return a pele Array which wraps the data in a numpy array
    
    Notes
    -----
    we must be careful that we only wrap the existing data
    """
    if not v.flags["FORC"]:
        raise ValueError("the numpy array is not c-contiguous.  copy it into a contiguous format before wrapping with pele::Array")
    return Array[size_t](<size_t *> v.data, v.size)

cdef inline Array[size_t] array_size_t_from_np(vin) except *:
    """return a pele Array which contains a copy of the data in a numpy array
    """
    cdef int i
    cdef np.ndarray[long, ndim=1] v = np.asarray(vin, dtype=long)
    cdef int N = v.size
    cdef Array[size_t] vnew = Array[size_t](N)
    for i in xrange(N):
        vnew[i] = v[i]
    return vnew

cdef inline np.ndarray[size_t, ndim=1] pele_array_to_np_size_t(Array[size_t] v):
    """copy the data in a pele::Array into a new numpy array
    """
    cdef int i
    cdef int N = v.size()
    cdef np.ndarray[size_t, ndim=1] vnew = np.zeros(N, dtype=size_t)
    for i in xrange(N):
        vnew[i] = v[i]
    return vnew
    
cdef inline Array[long] array_wrap_np_long(np.ndarray[long] v) except *:
    """return a pele Array which wraps the data in a numpy array
    
    Notes
    -----
    we must be careful that we only wrap the existing data
    """
    if not v.flags["FORC"]:
        raise ValueError("the numpy array is not c-contiguous.  copy it into a contiguous format before wrapping with pele::Array")
    return Array[long](<long*> v.data, v.size)
    
cdef inline np.ndarray[long, ndim=1] pele_array_to_np_long(Array[long] v):
    """copy the data in a pele::Array into a new numpy array
    """
    cdef int i
    cdef int N = v.size()
    cdef np.ndarray[long, ndim=1] vnew = np.zeros(N, dtype=long)
    for i in xrange(N):
        vnew[i] = v[i]
    return vnew
