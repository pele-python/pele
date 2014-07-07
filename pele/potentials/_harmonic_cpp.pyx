# distutils: language = c++
cimport pele.potentials._pele as _pele
import numpy as np
cimport numpy as np
#from cpython cimport bool 
from ctypes import c_size_t as size_t
from libcpp cimport bool as cbool

#===============================================================================
# THIS POTENTIAL NEEDS TO BE CLEANED UP
#===============================================================================

# use external c++ class
cdef extern from "pele/harmonic.h" namespace "pele":
    cdef cppclass cBaseHarmonic "pele::BaseHarmonic":
        cBaseHarmonic(_pele.Array[double] coords, double k, size_t ndim) except +
        void set_k(double) except +
        double get_k() except +
    cdef cppclass cHarmonic "pele::Harmonic":
        cHarmonic(_pele.Array[double] coords, double k, size_t ndim) except +
    cdef cppclass cHarmonicCOM "pele::HarmonicCOM":
        cHarmonicCOM(_pele.Array[double] coords, double k, size_t ndim) except +

cdef class Harmonic(_pele.BasePotential):
    """define the python interface to the c++ Harmonic potential implementation
    """
    cdef cBaseHarmonic* newptr
    cdef origin
    cdef k
    cdef pcom
    cdef pbdim
    
    def __cinit__(self, coords, pk, bdim=3, com=False):
        
        cdef np.ndarray[double, ndim=1] origin = coords
        cdef double k = pk
        cdef cbool pcom = com
        cdef int pbdim = bdim
        
        if pcom is True:
            self.thisptr = <_pele.cBasePotential*>new cHarmonicCOM(_pele.Array[double](<double*> origin.data, origin.size), k, pbdim)
        else:
            self.thisptr = <_pele.cBasePotential*>new cHarmonic(_pele.Array[double](<double*> origin.data, origin.size), k, pbdim)
            
        self.newptr = <cBaseHarmonic*> self.thisptr
        
        self.origin = origin
        self.k = k
        self.pcom = pcom
        self.pbdim = pbdim
        
    def set_k(self, newk):
        self.k = newk
        self.newptr.set_k(newk)
        
    def get_k(self):
        k = self.newptr.get_k()
        return k
    
    def __reduce__(self):
        return (Harmonic,(self.origin, self.k, self.pbdim, self.pcom))
        