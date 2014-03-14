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


#cdef extern from "_mypotential.c":
cdef extern void mypot_energy_gradient(double * coords, int * natoms, double * e,
                                       double * grad, double *eps, double *sig )


class MyPotIsoC(BasePotential):
    def __init__(self, natoms):
        self.natoms = natoms
    
    def getEnergyGradient(self, np.ndarray[double, ndim=1] coords):
        cdef np.ndarray[double, ndim=1] grad = np.zeros(coords.size)
#        cdef double * gdata = grad.data
#        cdef double * xdata = coords.data
        cdef double energy
        cdef double eps = 1.
        cdef double sig = 1.
        cdef int natoms = coords.size / 3
        mypot_energy_gradient(<double*> coords.data, <int *> &natoms, <double *> &energy, 
                    <double *>grad.data, <double*> &eps, <double*> &sig)
        return energy, grad
    
    def getEnergy(self, coords):
        e, g = self.getEnergyGradient(coords)
        return e
        
