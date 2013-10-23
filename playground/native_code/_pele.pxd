#cython: boundscheck=False
#cython: wraparound=False
##aaacython: noncheck=True

import numpy as np
cimport numpy as np 

cdef extern from "array.h" namespace "pele":
    cdef cppclass Array :
        Array() except +
        Array(double*, int n) except +
        size_t size()
        double *data()

cdef extern from "base_potential.h" namespace "pele":
    cdef cppclass  cBasePotential "pele::BasePotential":
        cBasePotential() except +
        double get_energy(Array &x) except *
        double get_energy_gradient(Array &x, Array &grad) except *
            
cdef extern from "potentialfunction.h" namespace "pele":
    cdef cppclass  cPotentialFunction "pele::PotentialFunction":
        cPotentialFunction(
            double (*energy)(Array x, void *userdata) except *,
            double (*energy_gradient)(Array x, Array grad, void *userdata) except *,
            void *userdata) except +
    
cdef class BasePotential:
    cdef cBasePotential *thisptr      # hold a C++ instance which we're wrapping
    
