# distutils: language = c++
import numpy as np

from ctypes import c_size_t as size_t
from libcpp cimport bool as cbool
cimport numpy as np

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np, array_wrap_np_size_t

#===============================================================================
# THIS POTENTIAL NEEDS TO BE CLEANED UP
#===============================================================================

# use external c++ class
cdef extern from "pele/pspin_spherical.h" namespace "pele":
    cdef cppclass cMeanFieldPSpinSpherical "pele::MeanFieldPSpinSpherical":
        cMeanFieldPSpinSpherical(_pele.Array[double] interaction, size_t nspins) except +

cdef class _Cdef_MeanFieldPSpinSpherical(_pele.BasePotential):
    """define the python interface to the c++ MeanFieldPSpinSpherical potential implementation
    """
    cdef cMeanFieldPSpinSpherical* newptr
    cdef interactions
    cdef int nspins
    
    def __cinit__(self, np.ndarray[double, ndim=1] interactions, int nspins):
        
        self.nspins = nspins
        self.interactions = interactions
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cMeanFieldPSpinSpherical(array_wrap_np(self.interactions), 
                                                                                                                       self.nspins) )
        self.newptr = <cMeanFieldPSpinSpherical*> self.thisptr.get()
    
    def __reduce__(self):
        return (MeanFieldPSpinSpherical,(self.interactions, self.nspins))

class MeanFieldPSpinSpherical(_Cdef_MeanFieldPSpinSpherical):
    """
    Python interface to _Cdef_MeanFieldPSpinSpherical
    """
        