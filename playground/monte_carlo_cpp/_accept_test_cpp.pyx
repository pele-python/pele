import numpy as np
cimport numpy as np
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport *
from _pele_mc import *

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef extern from "pele/accept_test.h" namespace "pele":
    cdef cppclass cppMetropolis "pele::Metropolis":
        cppMetropolis() except +

cdef class _Cdef_Metropolis(_Cdef_AcceptTest):
    """This class is the python interface for the c++ pele::Metropolis acceptance test class implementation
    """
    def __cinit__(self):
        self.thisptr = <cppAcceptTest*>new cppMetropolis()
        
        
class Metropolis(_Cdef_Metropolis):
    """This class is the python interface for the c++ Metropolis implementation.
    """