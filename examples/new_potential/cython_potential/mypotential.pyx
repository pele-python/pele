import numpy as np
cimport numpy as np

cimport cython

from pele.potentials import BasePotential

@cython.boundscheck(False)
@cython.cdivision(True)
cdef _getEnergyGradient(np.ndarray[double, ndim=2] coords, double sig, double eps):
    cdef double sig12 = sig**12
    cdef double sig24 = sig12**2
    cdef int natoms = len(coords[:,0])
    cdef double energy = 0.
    cdef np.ndarray[double, ndim=2] grad = np.zeros([natoms, 3], float)
    cdef np.ndarray[double, ndim=1] dr = np.zeros([3], float)
    cdef int i, j, k
    cdef double r2, ir2, g, ir12, ir24
    for i in xrange(natoms):
        for j in xrange(i):
            r2 = 0.
            for k in xrange(3):
                dr[k] = coords[i,k] - coords[j,k]
                r2 += dr[k] * dr[k]
            ir2 = 1. / r2
            ir12 = ir2 * ir2 * ir2
            ir12 *= ir12
            ir24 = ir12 * ir12
            energy += 4. * eps * ( sig24 * ir24 - sig12 * ir12 )
            g = 4. * eps * ir2 * ( -24. * sig24 * ir24 + 12. * sig12 * ir12)
            for k in xrange(3):
                grad[i,k] += g * dr[k]
                grad[j,k] -= g * dr[k]

    return energy, grad.reshape(-1)


class MyPotCython(BasePotential):
    def __init__(self, natoms, eps=1.0, sig=1.0):
        self.natoms = natoms
        self.eps = eps
        self.sig = sig

    def getEnergyGradient(self, coords):
        e, grad = _getEnergyGradient(coords.reshape([-1,3]), self.sig, self.eps)
        return e, grad
        
    def getEnergy(self, coords):
        e, grad = self.getEnergyGradient(coords)
        return e

