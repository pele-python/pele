cimport cython
import sys
from libcpp cimport bool as cbool
import numpy as np
cimport numpy as np
from pele.potentials import _pele
from pele.potentials cimport _pele
from pele.optimize cimport _pele_opt 
from _pele_mc cimport cppConfTest,_Cdef_ConfTest

cdef extern from "pele/conf_test.h" namespace "pele":
    cdef cppclass cppCheckSphericalContainer "pele::CheckSphericalContainer":
        cppCheckSphericalContainer(double) except+
    cdef cppclass cppCheckSameMinimum "pele::CheckSameMinimum":
        cppCheckSameMinimum(_pele_opt.cGradientOptimizer *, _pele.Array[double], _pele.Array[double] rattlers, 
                            double, double, double) except+

#===============================================================================
# Check spherical container
#===============================================================================

cdef class _Cdef_CheckSphericalContainer(_Cdef_ConfTest):
    """This class is the python interface for the c++ pele::CheckSphericalContainer configuration test class implementation
    """
    def __cinit__(self, radius):
        self.thisptr = <cppConfTest*>new cppCheckSphericalContainer(radius)
        
class CheckSphericalContainer(_Cdef_CheckSphericalContainer):
    """This class is the python interface for the c++ CheckSphericalContainer implementation.
    """

#===============================================================================
# Check same minimum
#===============================================================================

cdef class _Cdef_CheckSameMinimum(_Cdef_ConfTest):
    """This class is the python interface for the c++ pele::CheckSameMinimum configuration test class implementation
    """
    
    cdef _pele_opt.GradientOptimizer opt # this is stored so that the memory is not freed
    
    def __cinit__(self, optimizer, origin, rattlers, Eor, Etol, dtol):
        cdef np.ndarray[double, ndim=1] orginc = np.array(origin, dtype=float)
        cdef np.ndarray[double, ndim=1] rattlersc = np.array(rattlers, dtype=float)
        cdef _pele_opt.GradientOptimizer opt = optimizer
        self.thisptr = <cppConfTest*>new cppCheckSameMinimum(opt.thisptr, _pele.Array[double](<double*> orginc.data, orginc.size), 
                                                             _pele.Array[double](<double*> rattlersc.data, rattlersc.size), Eor, Etol, dtol)
        
class CheckSameMinimum(_Cdef_CheckSameMinimum):
    """This class is the python interface for the c++ CheckSameMinimum implementation.
    """