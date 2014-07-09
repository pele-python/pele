"""
# distutils: language = C++

basic potential interface stuff    
"""

import numpy as np
cimport numpy as np

# This is the base class for all potentials
cdef class BasePotential(object):
    """this class defines the python interface for c++ potentials 
    
    Notes
    -----
    for direct access to the underlying c++ potential use self.thisptr
    """
        
    def getEnergyGradient(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = x.copy()
        e = self.thisptr.get().get_energy_gradient(Array[double](<double*> x.data, x.size),
                                             Array[double](<double*> grad.data, grad.size))
        return e, grad
    
    def getEnergy(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        return self.thisptr.get().get_energy(Array[double](<double*> x.data, x.size))
    
    def getGradient(self, np.ndarray[double, ndim=1] x not None):
        e, grad = self.getEnergyGradient(x)
        return grad
    
    def getEnergyGradientHessian(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x.size)
        cdef np.ndarray[double, ndim=1] hess = np.zeros(x.size**2)
        e = self.thisptr.get().get_energy_gradient_hessian(Array[double](<double*> x.data, x.size),
                                             Array[double](<double*> grad.data, grad.size),
                                             Array[double](<double*> hess.data, hess.size),
                                             )
        return e, grad, hess.reshape([x.size, x.size])
    
    def getHessian(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] hess = np.zeros(x.size**2)
        self.thisptr.get().get_hessian(Array[double](<double*> x.data, x.size),
                                 Array[double](<double*> hess.data, hess.size),
                                 )
        return np.reshape(hess, [x.size, x.size])
    
    def NumericalDerivative(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = np.zeros([x.size])
        self.thisptr.get().numerical_gradient(Array[double](<double*> x.data, x.size),
                                       Array[double](<double*> grad.data, grad.size),
                                       eps
                                       )
        return grad
                
    def NumericalHessian(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] hess = np.zeros([x.size**2])
        self.thisptr.get().numerical_hessian(Array[double](<double*> x.data, x.size),
                                       Array[double](<double*> hess.data, hess.size),
                                       eps
                                       )
#        newhess = hess;
        return np.reshape(hess, [x.size, x.size])

# This is a little test function to benchmark potential evaluation in a loop
# in native code    
#def call_pot(Potential pot, np.ndarray[double, ndim=1, mode="c"] x not None,
#              np.ndarray[double, ndim=1, mode="c"] grad not None,
#              int N):
#    _call_pot(pot.thisptr, Array[double](<double*> x.data, x.size),
#             Array[double](<double*> grad.data, grad.size), N)
