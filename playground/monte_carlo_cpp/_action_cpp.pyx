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

cdef extern from "pele/actions.h" namespace "pele":
    cdef cppclass cppRecordEnergyHistogram "pele::RecordEnergyHistogram":
        cppRecordEnergyHistogram(double, double, double) except +
        void get_histogram(_pele.Array[size_t]&) except +

cdef class _Cdef_RecordEnergyHistogram(_Cdef_Action):
    """This class is the python interface for the c++ pele::RecordEnergyHistogram acceptance test class implementation
    """
    def __cinit__(self, min, max, bin):
        self.thisptr = <cppAction*>new cppRecordEnergyHistogram(min, max, bin)
    
    def get_histogram(self, hist):
        cdef np.ndarray[double, ndim=1] histc = np.array(hist, dtype=float)
        self.thisptr.get_histogram(histc)
        
class RecordEnergyHistogram(_Cdef_RecordEnergyHistogram):
    """This class is the python interface for the c++ RecordEnergyHistogram implementation.
    """