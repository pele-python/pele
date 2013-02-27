#cython: boundscheck=False
#cython: wraparound=False
##aaacython: noncheck=True

import numpy as np
cimport numpy as np 

cdef extern from "potential.h" namespace "pygmin":
    cdef cppclass Array :
        Array(double*, int n) except +

    cdef cppclass  cPotential "pygmin::Potential":
        cPotential() except +
        double get_energy(Array &x)
        double get_energy_gradient(Array &x, Array &grad)

    cdef cppclass  cPotentialFunction "pygmin::PotentialFunction":
        cPotentialFunction(
        	double (*energy)(double *x, int n, void *userdata),
            double (*energy_gradient)(double *x, double *grad, int n, void *userdata),
            void *userdata) except +
            
    cdef void _call_pot "call_pot" (cPotential *pot, Array &x, Array &grad, int n)
    
cdef class Potential:
    cdef cPotential *thisptr      # hold a C++ instance which we're wrapping
    