"""
# distutils: language = C++

basic potential interface stuff    
"""
import numpy as np
cimport numpy as np
from ctypes import c_size_t as size_t


cdef class BasePotential(object):
    """this class defines the python interface for c++ potentials 
    
    Notes
    -----
    for direct access to the underlying c++ potential use self.thisptr
    """
        
    def getEnergyGradient(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x.size)
        e = self.thisptr.get().get_energy_gradient(array_wrap_np(x),
                                                   array_wrap_np(grad))
        return e, grad
    
    def getEnergy(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        return self.thisptr.get().get_energy(array_wrap_np(x))
    
    def getGradient(self, np.ndarray[double, ndim=1] x not None):
        e, grad = self.getEnergyGradient(x)
        return grad
    
    def getEnergyGradientHessian(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x.size)
        cdef np.ndarray[double, ndim=1] hess = np.zeros(x.size**2)
        e = self.thisptr.get().get_energy_gradient_hessian(array_wrap_np(x),
                                                           array_wrap_np(grad),
                                                           array_wrap_np(hess))
        return e, grad, hess.reshape([x.size, x.size])
    
    def getHessian(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] hess = np.zeros(x.size**2)
        self.thisptr.get().get_hessian(array_wrap_np(x), array_wrap_np(hess))
        return np.reshape(hess, [x.size, x.size])
    
    def NumericalDerivative(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = np.zeros([x.size])
        self.thisptr.get().numerical_gradient(array_wrap_np(x),
                                              array_wrap_np(grad), 
                                              eps)
        return grad
                
    def NumericalHessian(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] hess = np.zeros([x.size**2])
        self.thisptr.get().numerical_hessian(array_wrap_np(x),
                                             array_wrap_np(hess),
                                             eps)
        return np.reshape(hess, [x.size, x.size])
