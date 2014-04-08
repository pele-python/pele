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
    
    Notes
    -----
    for direct access to the underlying c++ potential use self.thisptr
    """
    def __cinit__(self, *args, **kwargs):
        # store an instance to the current c++ class, will be used in every call
        self.thisptr = NULL#<cBasePotential*>new cBasePotential()
    
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    
#    def getEnergyGradientInplace(self,
#                        np.ndarray[double, ndim=1] x not None,
#                        np.ndarray[double, ndim=1] grad not None):
#        # redirect the call to the c++ class
#        e = self.thisptr.get_energy_gradient(Array[double](<double*> x.data, x.size),
#                                             Array[double](<double*> grad.data, grad.size))
#        return e
        
    def getEnergyGradient(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = x.copy()
        e = self.thisptr.get_energy_gradient(Array[double](<double*> x.data, x.size),
                                             Array[double](<double*> grad.data, grad.size))
        return e, grad
    
    def getEnergy(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        return self.thisptr.get_energy(Array[double](<double*> x.data, x.size))
    
    def getGradient(self, np.ndarray[double, ndim=1] x not None):
        e, grad = self.getEnergyGradient(x)
        return grad
    
    def getEnergyGradientHessian(self, np.ndarray[double, ndim=1] x not None):
        e, grad = self.getEnergyGradient(x)
        hess = self.Hessian(x) #this should return the numerical hessian if hessian is not implemented
        return e, grad, hess
    
    def getHessian(self, np.ndarray[double, ndim=1] x not None):
        e, grad, hess = self.getEnergyGradientHessian(x)
        return hess
    
    def NumericalDerivative(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = np.zeros([x.size])
        self.thisptr.numerical_gradient(Array[double](<double*> x.data, x.size),
                                       Array[double](<double*> grad.data, grad.size),
                                       eps
                                       )
        return grad
                
    def Hessian(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] hess = np.zeros([x.size**2])
        self.thisptr.get_hessian(Array[double](<double*> x.data, x.size),
                                       Array[double](<double*> hess.data, hess.size))
#        newhess = hess;
        return np.reshape(hess, [x.size, x.size])    
        
    def NumericalHessian(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] hess_num = np.zeros([x.size**2])
        self.thisptr.numerical_hessian(Array[double](<double*> x.data, x.size),
                                       Array[double](<double*> hess_num.data, hess_num.size),
                                       eps
                                       )
#        newhess = hess;
        return np.reshape(hess_num, [x.size, x.size])

# This is a little test function to benchmark potential evaluation in a loop
# in native code    
#def call_pot(Potential pot, np.ndarray[double, ndim=1, mode="c"] x not None,
#              np.ndarray[double, ndim=1, mode="c"] grad not None,
#              int N):
#    _call_pot(pot.thisptr, Array[double](<double*> x.data, x.size),
#             Array[double](<double*> grad.data, grad.size), N)
