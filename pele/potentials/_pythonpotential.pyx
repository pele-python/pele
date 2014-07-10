"""
# distutils: language = C++
"""
import numpy as np

cimport numpy as np
cimport numpy as cnp
cimport cython
from cpython.ref cimport PyObject

from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr

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
    
cdef class CppPotentialWrapperBase(_pele.BasePotential):
    def __cinit__(self, *args, **kwargs):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cPythonPotential(
                                           <PyObject*>self) )

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
        
