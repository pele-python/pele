import numpy as np
cimport numpy as np
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport cppAction,_Cdef_Action

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef extern from "pele/actions.h" namespace "pele":
    cdef cppclass cppRecordEnergyHistogram "pele::RecordEnergyHistogram":
        cppRecordEnergyHistogram(double, double, double) except +
        void get_histogram(_pele.Array[double]&) except +
        size_t get_histogram_size() except +
        void print_histogram(size_t) except +
        
cdef class _Cdef_RecordEnergyHistogram(_Cdef_Action):
    """This class is the python interface for the c++ pele::RecordEnergyHistogram acceptance test class implementation
    """
    def __cinit__(self, min, max, bin):
        self.thisptr = <cppAction*>new cppRecordEnergyHistogram(min, max, bin)
    
    def get_histogram(self, np.ndarray[double, ndim=1] hist):
        cdef cppRecordEnergyHistogram* newptr = <cppRecordEnergyHistogram*> self.thisptr
        newptr.get_histogram(_pele.Array[double](<double*> hist.data, hist.size))
        
    def get_histogram_size(self):
        cdef cppRecordEnergyHistogram* newptr2 = <cppRecordEnergyHistogram*> self.thisptr
        size = newptr2.get_histogram_size()
        return size
    
    def print_histogram(self, ntot):
        cdef cppRecordEnergyHistogram* newptr3 = <cppRecordEnergyHistogram*> self.thisptr
        newptr3.print_histogram(ntot)
        
class RecordEnergyHistogram(_Cdef_RecordEnergyHistogram):
    """This class is the python interface for the c++ RecordEnergyHistogram implementation.
    """