from playground.native_code cimport _pele
import numpy as np
cimport numpy as np

# use external c++ class
cdef extern from "lj.h" namespace "pele":
    cdef cppclass  cLJ "pele::LJ":
        cLJ(double C6, double C12) except +
    cdef cppclass  cLJ_Ilist "pele::LJ_interaction_list":
        cLJ_Ilist(_pele.Array[np.int64_t] & ilist, double C6, double C12) except +

cdef class LJ(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    def __cinit__(self, eps=1.0, sigma=1.0):
        self.thisptr = <_pele.cBasePotential*>new cLJ(4.*eps*sigma**6, 4.*eps*sigma**12)


cdef class LJInteractionList(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    def __cinit__(self, np.ndarray[np.int64_t, ndim=1] ilist, eps=1.0, sigma=1.0):
        self.thisptr = <_pele.cBasePotential*>new cLJ_Ilist( _pele.Array[np.int64_t](<np.int64_t*> ilist.data, <int> ilist.size), 4.*eps*sigma**6, 4.*eps*sigma**12)
