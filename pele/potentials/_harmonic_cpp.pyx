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

cdef class Harmonic(_pele.BasePotential):
    """define the python interface to the c++ Harmonic potential implementation
    """
    #thickness of the wca shell is sca * R where R is the hard core radius of the sphere
    def __cinit__(self, np.ndarray[double, ndim=1] coords, double k, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
#        if boxl is not None:
#            boxvec = [boxl] * 3
#        cdef np.ndarray[double, ndim=1] bv
                
        #if boxvec is None:
        #self.periodic = False
        self.thisptr = <_pele.cBasePotential*>new cHarmonic(_pele.Array[double](<double*> coords.data, coords.size),k)
        #else:
        #raise NameError('periodic potentialharmonic potential')
    
    def set_k(self, newk):
        cdef cHarmonic* newptr = <cHarmonic*> self.thisptr
        newptr.set_k(newk)
  
