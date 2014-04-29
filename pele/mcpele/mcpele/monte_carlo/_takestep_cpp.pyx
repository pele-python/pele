# distutils: language = c++
cimport cython
import sys
from _pele_mc cimport cppTakeStep,_Cdef_TakeStep

#===============================================================================
# RandomCoordsDisplacement
#===============================================================================

cdef extern from "pele/takestep.h" namespace "pele":
    cdef cppclass cppRandomCoordsDisplacement "pele::RandomCoordsDisplacement":
        cppRandomCoordsDisplacement(size_t, size_t) except +

cdef class _Cdef_RandomCoordsDisplacement(_Cdef_TakeStep):
    """This class is the python interface for the c++ pele::RandomCoordsDisplacement take step class implementation
    """
    def __cinit__(self, ndim, rseed):
        self.thisptr = <cppTakeStep*>new cppRandomCoordsDisplacement(ndim, rseed)

class RandomCoordsDisplacement(_Cdef_RandomCoordsDisplacement):
    """This class is the python interface for the c++ RandomCoordsDisplacement implementation.
    """