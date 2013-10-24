#cython: boundscheck=False
#cython: wraparound=False
##aaacython: noncheck=True

import numpy as np
cimport numpy as np 

cdef extern from "array.h" namespace "pele":
    cdef cppclass Array[dtype] :
        Array() except +
        Array(dtype*, int n) except +
        size_t size()
        dtype *data()

cdef extern from "base_potential.h" namespace "pele":
    cdef cppclass  cBasePotential "pele::BasePotential":
        cBasePotential() except +
        double get_energy(Array[double] &x) except *
        double get_energy_gradient(Array[double] &x, Array[double] &grad) except *
            
cdef extern from "potentialfunction.h" namespace "pele":
    cdef cppclass  cPotentialFunction "pele::PotentialFunction":
        cPotentialFunction(
            double (*energy)(Array[double] x, void *userdata) except *,
            double (*energy_gradient)(Array[double] x, Array[double] grad, void *userdata) except *,
            void *userdata) except +
    
cdef class BasePotential:
    cdef cBasePotential *thisptr      # hold a C++ instance which we're wrapping
    
