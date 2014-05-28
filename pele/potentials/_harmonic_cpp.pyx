# distutils: language = c++
cimport pele.potentials._pele as _pele
import numpy as np
cimport numpy as np
from cpython cimport bool
from ctypes import c_size_t as size_t

#===============================================================================
# THIS POTENTIAL NEEDS TO BE CLEANED UP
#===============================================================================

# use external c++ class
cdef extern from "pele/harmonic.h" namespace "pele":
    cdef cppclass  cBaseHarmonic "pele::BaseHarmonic":
        cBaseHarmonic(_pele.Array[double] coords, double k, size_t ndim) except +
        void set_k(double) except +
        double get_k() except +
    cdef cppclass  cHarmonic "pele::Harmonic":
        cHarmonic(_pele.Array[double] coords, double k, size_t ndim) except +
    cdef cppclass  cHarmonicCOM "pele::HarmonicCOM":
        cHarmonicCOM(_pele.Array[double] coords, double k, size_t ndim) except +

cdef class Harmonic(_pele.BasePotential):
    """define the python interface to the c++ Harmonic potential implementation
    """
    cdef cBaseHarmonic* newptr
    cdef origin
    cdef k
    cdef pcom
    cdef pndim
    
    def __cinit__(self, np.ndarray[double, ndim=1] coords, double pk, com, ndim):
        
        cdef np.ndarray[double, ndim=1] origin = coords
        cdef double k = pk
        cdef bool pcom = com
        cdef int pndim = ndim
        
        if com is True:
            self.thisptr = <_pele.cBasePotential*>new cHarmonicCOM(_pele.Array[double](<double*> origin.data, origin.size), k, pndim)
        else:
            self.thisptr = <_pele.cBasePotential*>new cHarmonic(_pele.Array[double](<double*> origin.data, origin.size), k, pndim)
            
        self.newptr = <cBaseHarmonic*> self.thisptr
        
        self.origin = origin
        self.k = k
        self.pcom = pcom
        self.pndim = pndim
        
    def set_k(self, newk):
        self.k = newk
        self.newptr.set_k(newk)
        
    def get_k(self):
        k = self.newptr.get_k()
        return k
    
    def __reduce__(self):
        return (Harmonic,(self.origin, self.k, self.pcom, self.pndim))
        