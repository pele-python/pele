# distutils: language = c++
import numpy as np

from ctypes import c_size_t as size_t
from libcpp cimport bool as cbool
cimport numpy as np

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np, array_wrap_np_size_t

# cython has no support for integer template argument.  This is a hack to get around it
# https://groups.google.com/forum/#!topic/cython-users/xAZxdCFw6Xs
# Basically you fool cython into thinking INT2 is the type integer,
# but in the generated c++ code you use 2 instead.
# The cython code MyClass[INT2] will create c++ code MyClass<2>.
cdef extern from *:
    ctypedef int INT2 "2"    # a fake type
    ctypedef int INT3 "3"    # a fake type
    ctypedef int INT4 "4"    # a fake type
    ctypedef int INT5 "5"    # a fake type
    
#===============================================================================
# THIS POTENTIAL NEEDS TO BE CLEANED UP
#===============================================================================

# use external c++ class
cdef extern from "pele/pspin_spherical.h" namespace "pele":
    cdef cppclass cMeanFieldPSpinSpherical "pele::MeanFieldPSpinSpherical"[p]:
        cMeanFieldPSpinSpherical(_pele.Array[double] interaction, size_t nspins, double tol) except +

cdef class _Cdef_MeanFieldPSpinSpherical(_pele.BasePotential):
    """define the python interface to the c++ MeanFieldPSpinSpherical potential implementation
    """
    cdef interactions
    cdef int nspins
    cdef int p
    cdef double tol

    def __cinit__(self, np.ndarray[double, ndim=1] interactions, int nspins, int p, double tol=1e-7):
        
        self.nspins = nspins
        self.p = p
        self.interactions = interactions
        self.tol = tol
        assert interactions.size == np.power(nspins, p), "MeanFieldPSpinSpherical: interactions is of the wrong size"
        
        if self.p == 2:
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cMeanFieldPSpinSpherical[INT2](array_wrap_np(self.interactions), 
                                                                                                                 self.nspins, self.tol) )
        elif self.p == 3:
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cMeanFieldPSpinSpherical[INT3](array_wrap_np(self.interactions), 
                                                                                                                 self.nspins, self.tol) )
        elif self.p == 4:
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cMeanFieldPSpinSpherical[INT4](array_wrap_np(self.interactions), 
                                                                                                                 self.nspins, self.tol) )
        elif self.p == 5:
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cMeanFieldPSpinSpherical[INT5](array_wrap_np(self.interactions), 
                                                                                                                 self.nspins, self.tol) )
        else:
            raise Exception("MeanFieldPSpinSpherical: illegal p")
            
    def __reduce__(self):
        return (MeanFieldPSpinSpherical,(self.interactions, self.nspins, self.p, self.tol))

class MeanFieldPSpinSpherical(_Cdef_MeanFieldPSpinSpherical):
    """
    Python interface to _Cdef_MeanFieldPSpinSpherical
    """