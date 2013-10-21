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

cdef extern from "potential.h" namespace "pele":
    cdef cppclass  cPotential "pele::Potential":
        cPotential() except +
        double get_energy(Array &x)
        double get_energy_gradient(Array &x, Array &grad)
            
cdef extern from "potentialfunction.h" namespace "pele":
    cdef cppclass  cPotentialFunction "pele::PotentialFunction":
        cPotentialFunction(
        	double (*energy)(Array x, void *userdata),
            double (*energy_gradient)(Array x, Array grad, void *userdata),
            void *userdata) except +
    
cdef class Potential:
    cdef cPotential *thisptr      # hold a C++ instance which we're wrapping
    
