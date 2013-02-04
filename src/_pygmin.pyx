#cython: boundscheck=False
#cython: wraparound=False
##aaacython: noncheck=True

import numpy as np
cimport numpy as np

cdef class Potential:
    def __cinit__(self):
        self.thisptr = <cPotential*>new cPotential()
    
    def __dealloc__(self):
        del self.thisptr
    
    def get_energy_gradient_inplace(self,
                        np.ndarray[double, ndim=1] x not None,
                        np.ndarray[double, ndim=1] grad not None):
        e = self.thisptr.get_energy_gradient(Array(<double*> x.data, x.size),
                                             Array(<double*> grad.data, grad.size))
        return e
        
    def get_energy_gradient(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] grad = x.copy()
        e = self.thisptr.get_energy_gradient(Array(<double*> x.data, x.size),
                                             Array(<double*> grad.data, grad.size))
        return e, grad
            
def call_pot(Potential pot, np.ndarray[double, ndim=1, mode="c"] x not None,
              np.ndarray[double, ndim=1, mode="c"] grad not None,
              int N):
    _call_pot(pot.thisptr, Array(<double*> x.data, x.size),
             Array(<double*> grad.data, grad.size), N)
