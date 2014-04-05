import numpy as np
cimport numpy as np
from pele.potentials import _pele
from pele.potentials cimport _pele
cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport cppAction,_Cdef_Action

#===============================================================================
# Adjust step size
#===============================================================================

cdef extern from "pele/actions.h" namespace "pele": 
    cdef cppclass cppAdjustStep "pele::AdjustStep":
        cppAdjustStep(double, double, size_t, size_t) except +
        
cdef class _Cdef_AdjustStep(_Cdef_Action):
    """This class is the python interface for the c++ pele::AdjustStep action class implementation
    """
    def __cinit__(self, target, factor, niter, navg):
        self.thisptr = <cppAction*>new cppAdjustStep(target, factor, niter, navg)
        
class AdjustStep(_Cdef_AdjustStep):
    """This class is the python interface for the c++ AdjustStep implementation.
    """

#===============================================================================
# Record Energy Histogram
#===============================================================================        

cdef extern from "pele/actions.h" namespace "pele":
    cdef cppclass cppRecordEnergyHistogram "pele::RecordEnergyHistogram":
        cppRecordEnergyHistogram(double, double, double, size_t) except +
        _pele.Array[double] get_histogram() except +
        void print_terminal(size_t) except +
        double get_max() except +
        double get_min() except +
        
cdef class _Cdef_RecordEnergyHistogram(_Cdef_Action):
    """This class is the python interface for the c++ pele::RecordEnergyHistogram acceptance test class implementation
    """
    def __cinit__(self, min, max, bin, eqsteps):
        self.thisptr = <cppAction*>new cppRecordEnergyHistogram(min, max, bin, eqsteps)
    
    @cython.boundscheck(False)
    def get_histogram(self):
        """return a histogram array"""
        cdef cppRecordEnergyHistogram* newptr = <cppRecordEnergyHistogram*> self.thisptr
        cdef _pele.Array[double] histi = newptr.get_histogram()
        cdef double *histdata = histi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] hist = np.zeros(histi.size())
        cdef size_t i
        for i in xrange(histi.size()):
            hist[i] = histdata[i]
              
        return hist
        
    def print_terminal(self, ntot):
        cdef cppRecordEnergyHistogram* newptr2 = <cppRecordEnergyHistogram*> self.thisptr
        newptr2.print_terminal(ntot)
    
    def get_bounds_val(self):
        cdef cppRecordEnergyHistogram* newptr3 = <cppRecordEnergyHistogram*> self.thisptr
        Emin = newptr3.get_min()
        Emax = newptr3.get_max()
        return Emin, Emax
        
class RecordEnergyHistogram(_Cdef_RecordEnergyHistogram):
    """This class is the python interface for the c++ RecordEnergyHistogram implementation.
    """
    
#===============================================================================
# RecordDisp2Histogram
#===============================================================================
#derives from record energy histogram
cdef extern from "pele/actions.h" namespace "pele":
    cdef cppclass cppRecordDisp2Histogram "pele::RecordDisp2Histogram":
        cppRecordDisp2Histogram(_pele.Array[double],_pele.Array[double],double, double, double, size_t) except +
        _pele.Array[double] get_histogram() except +
        void print_terminal(size_t) except +
        double get_max() except +
        double get_min() except +
    cdef cppclass cppRecordDisp2HistogramPeriodic "pele::RecordDisp2HistogramPeriodic":
        cppRecordDisp2HistogramPeriodic(_pele.Array[double],_pele.Array[double],double, double, double, size_t, double * boxvec) except +
        _pele.Array[double] get_histogram() except +
        void print_terminal(size_t) except +
        double get_max() except +
        double get_min() except +

cdef class _Cdef_RecordDisp2Histogram(_Cdef_Action):
    """This class is the python interface for the c++ pele::RecordDisp2Histogram acceptance test class implementation
    """
    cpdef cbool periodic 
    def __cinit__(self, origin, rattlers, min, max, bin, eqsteps, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * 3
        cdef np.ndarray[double, ndim=1] bv
        
        cdef np.ndarray[double, ndim=1] orginc = np.array(origin, dtype=float)
        cdef np.ndarray[double, ndim=1] rattlersc = np.array(rattlers, dtype=float)
        
        if boxvec is None:
            self.periodic = False
            self.thisptr = <cppAction*>new cppRecordDisp2Histogram(_pele.Array[double](<double*> orginc.data, orginc.size),
                                                               _pele.Array[double](<double*> rattlersc.data, rattlersc.size),
                                                               min, max, bin, eqsteps)
        else:
            self.periodic = True
            bv = np.array(boxvec, dtype=float)
            self.thisptr = <cppAction*>new cppRecordDisp2HistogramPeriodic(_pele.Array[double](<double*> orginc.data, orginc.size),
                                                               _pele.Array[double](<double*> rattlersc.data, rattlersc.size),
                                                               min, max, bin, eqsteps, <double*> bv.data)
    @cython.boundscheck(False)
    def get_histogram(self):
        """return a histogram array"""
        cdef cppRecordEnergyHistogram* newptr = <cppRecordEnergyHistogram*> self.thisptr
        cdef _pele.Array[double] histi = newptr.get_histogram()
        cdef double *histdata = histi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] hist = np.zeros(histi.size())
        cdef size_t i
        for i in xrange(histi.size()):
            hist[i] = histdata[i]
              
        return hist
        
    def print_terminal(self, ntot):
        cdef cppRecordEnergyHistogram* newptr2 = <cppRecordEnergyHistogram*> self.thisptr
        newptr2.print_terminal(ntot)
    
    def get_bounds_val(self):
        cdef cppRecordEnergyHistogram* newptr3 = <cppRecordEnergyHistogram*> self.thisptr
        dmin = newptr3.get_min()
        dmax = newptr3.get_max()
        return dmin, dmax
        
class RecordDisp2Histogram(_Cdef_RecordDisp2Histogram):
    """This class is the python interface for the c++ RecordDisp2Histogram implementation.
    """