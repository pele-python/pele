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
    cdef cppclass  cHarmonic "pele::Harmonic":
        cHarmonic(_pele.Array[double] coords, double k) except +
        void set_k(double) except +
        double get_k() except +
    cdef cppclass  cHarmonicPeriodic "pele::HarmonicPeriodic":
        cHarmonicPeriodic(_pele.Array[double] coords, double k, double * boxvec) except +
        void set_k(double) except +
        double get_k() except +

cdef class Harmonic(_pele.BasePotential):
    """define the python interface to the c++ Harmonic potential implementation
    """
    cpdef bool periodic 
    cdef cHarmonic* newptr
    def __cinit__(self, np.ndarray[double, ndim=1] coords, double k, boxvec=None, boxl=None):
        self.newptr = <cHarmonic*> self.thisptr
        
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * 3
        cdef np.ndarray[double, ndim=1] bv
                
        if boxvec is None:
            self.periodic = False
            self.thisptr = <_pele.cBasePotential*>new cHarmonic(_pele.Array[double](<double*> coords.data, coords.size),k)
        else:
            self.periodic = True
            bv = np.array(boxvec, dtype=float)
            self.thisptr = <_pele.cBasePotential*>new cHarmonicPeriodic(_pele.Array[double](<double*> coords.data, coords.size), k, <double*> bv.data)
    
    def set_k(self, newk):
        self.newptr.set_k(newk)
        
    def get_k(self):
        k = self.newptr.get_k()
        return k
  
