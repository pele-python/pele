import numpy as np
cimport numpy as np
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport cppConfTest,_Cdef_ConfTest

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef extern from "pele/conf_test.h" namespace "pele":
    cdef cppclass cppCheckSphericalContainer "pele::CheckSphericalContainer":
        cppCheckSphericalContainer(double) except +

cdef class _Cdef_CheckSphericalContainer(_Cdef_ConfTest):
    """This class is the python interface for the c++ pele::CheckSphericalContainer configuration test class implementation
    """
    def __cinit__(self, radius):
        self.thisptr = <cppConfTest*>new cppCheckSphericalContainer(radius)
        
        
class CheckSphericalContainer(_Cdef_CheckSphericalContainer):
    """This class is the python interface for the c++ CheckSphericalContainer implementation.
    """