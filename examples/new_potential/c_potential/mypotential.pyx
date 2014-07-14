"""
This file must first be `cythonized` to a c file mypotential.c

    cython mypotential.pyx

We now use setup.py to create the shared object library

    python setup.py build_ext -i
"""
cimport cython
import numpy as np
cimport numpy as np
from pele.potentials import BasePotential


cdef extern from "_mypotential.h":
    double mypotential(double *x, int N, double *grad, double eps, double sig)


class MyPotC(BasePotential):
    def __init__(self, natoms):
        self.natoms = natoms
    
    def getEnergyGradient(self, np.ndarray[double, ndim=1] coords):
        cdef np.ndarray[double, ndim=1] grad = np.zeros(coords.size)
#        cdef double * gdata = grad.data
#        cdef double * xdata = coords.data
        cdef double energy = mypotential(<double*> coords.data, <int> coords.size, <double *>grad.data, 1., 1.)
        return energy, grad
    
    def getEnergy(self, coords):
        e, g = self.getEnergyGradient(coords)
        return e
        
