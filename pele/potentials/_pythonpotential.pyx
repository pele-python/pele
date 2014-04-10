"""
# distutils: language = C++
"""

from pele.potentials cimport _pele
import numpy as np
cimport numpy as np
cimport numpy as cnp
cimport cython
from cpython.ref cimport PyObject

# import the numpy array api.  this is necessary because some files
# used by this module use the numpy c-api.
# see the miscillaneous section and "importing the api"
# http://docs.scipy.org/doc/numpy/reference/c-api.array.html
# also see http://mail.scipy.org/pipermail/numpy-discussion/2011-December/059612.html
# question: should I define PY_ARRAY_UNIQUE_SYMBOL?  And how can I?
cnp.import_array()

cdef extern from "pele/python_potential_wrapper.h" namespace "pele":
    cdef cppclass  cPythonPotential "pele::PythonPotential":
        cPythonPotential(PyObject *potential) except +
    
# @cython.boundscheck(False)
# cdef double _python_grad(_pele.Array[double] x, _pele.Array[double] grad, void *userdata) except *:
#     cdef double *xdata = x.data()
#     cdef double *gdata = grad.data()
#     cdef np.ndarray[double, ndim=1, mode="c"] px = np.zeros(x.size())
#     cdef size_t i
#     for i in xrange(x.size()):
#         px[i] = xdata[i]
# 
# #    print "casting as pythonpotential"
#     pot = <PythonPotential>(userdata)
# #    print "computing energy"
#     
#     cdef np.ndarray[double, ndim=1] g
#     e, g = pot.getEnergyGradient(px)
# 
#     
#     for i in xrange(x.size()):
#         gdata[i] = g[i]
#     
# #    print "returning"
#     return e
# 
# # energy callback not yet implemented
# @cython.boundscheck(False)
# cdef double _python_energy(_pele.Array[double] x, void *userdata) except *:
#     cdef double *xdata = x.data()
#     cdef np.ndarray[double, ndim=1, mode="c"] px = np.zeros(x.size())
#     cdef size_t i
#     for i in xrange(x.size()):
#         px[i] = xdata[i]
# 
#     pot = <PythonPotential>(userdata)
#     return pot.getEnergy(px)

# define the potential class
# cdef class PythonPotential(_pele.BasePotential):   
#     def __cinit__(self, *args, **kwargs):
#         self.thisptr = <_pele.cBasePotential*>new _pele.cPotentialFunction(
#                                            &_python_energy,
#                                            &_python_grad,
#                                            <void*>self)
#     
#     def __dealloc__(self):
#         if self.thisptr != NULL:
#             del self.thisptr
#             self.thisptr = NULL
#     
#     def getEnergy(self, x):
#         raise NotImplementedError        
# 
#     def getEnergyGradient(self, x):
#         raise NotImplementedError
# 
# class CppPotentialWrapperOld(PythonPotential):
#     """wrap a python potential to be used in c++"""
#     def __init__(self, pot):
#         self.pot = pot
#         self.getEnergy = pot.getEnergy
#         self.getEnergyGradient = pot.getEnergyGradient

cdef class CppPotentialWrapperBase(_pele.BasePotential):
#     cdef potential 
    def __cinit__(self, *args, **kwargs):
        self.thisptr = <_pele.cBasePotential*>new cPythonPotential(
                                           <PyObject*>self)
#         self.potential = potential
#         self.getEnergy = potential.getEnergy
#         self.getEnergyGradient = potential.getEnergyGradient

class CppPotentialWrapper(CppPotentialWrapperBase):
    """wrap a python potential to be used in c++"""
    def __init__(self, pot):
        self.pot = pot
        self.getEnergy = pot.getEnergy
        self.getEnergyGradient = pot.getEnergyGradient
#         self.NumericalGradient = pot.NumericalGradient


class _TestingCppPotentialWrapper(CppPotentialWrapper):
    """testing potential which provides direct access to c++ wrapper"""
    def getEnergy(self, x):
        print "going through cpp"
        return _pele.BasePotential.getEnergy(self, x)

    def cpp_get_energy(self, x):
        return _pele.BasePotential.getEnergy(self, x)

    def cpp_get_energy_gradient(self, x):
        return _pele.BasePotential.getEnergyGradient(self, x)

    def getEnergyGradient(self, x):
        return _pele.BasePotential.getEnergyGradient(self, x)
        
