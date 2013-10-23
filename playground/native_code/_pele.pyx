#cython: boundscheck=False
#cython: wraparound=False
##aaacython: noncheck=True
#
#
# basic potential interface stuff    
#
import numpy as np
cimport numpy as np

# This is the base class for all potentials
cdef class BasePotential(object):
    """this class defines the python interface for c++ potentials 
    """
    def __cinit__(self):
        # store an instance to the current c++ class, will be used in every call
        self.thisptr = <cBasePotential*>new cBasePotential()
    
    def __dealloc__(self):
        del self.thisptr
    
    def getEnergyGradientInplace(self,
                        np.ndarray[double, ndim=1] x not None,
                        np.ndarray[double, ndim=1] grad not None):
        # redirect the call to the c++ class
        e = self.thisptr.get_energy_gradient(Array(<double*> x.data, x.size),
                                             Array(<double*> grad.data, grad.size))
        return e
        
    def getEnergyGradient(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = x.copy()
        e = self.thisptr.get_energy_gradient(Array(<double*> x.data, x.size),
                                             Array(<double*> grad.data, grad.size))
        return e, grad
    
    def getEnergy(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        return self.thisptr.get_energy(Array(<double*> x.data, x.size))
                
# This is a little test function to benchmark potential evaluation in a loop
# in native code    
#def call_pot(Potential pot, np.ndarray[double, ndim=1, mode="c"] x not None,
#              np.ndarray[double, ndim=1, mode="c"] grad not None,
#              int N):
#    _call_pot(pot.thisptr, Array(<double*> x.data, x.size),
#             Array(<double*> grad.data, grad.size), N)
