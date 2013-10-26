from playground.native_code cimport _pele
import numpy as np
cimport numpy as np
from cpython cimport bool

# use external c++ class
cdef extern from "lj.h" namespace "pele":
    cdef cppclass  cLJ "pele::LJ":
        cLJ(double C6, double C12) except +
    cdef cppclass  cLJPeriodic "pele::LJPeriodic":
        cLJPeriodic(double C6, double C12, double * boxvec) except +
    cdef cppclass  cLJCut "pele::LJCut":
        cLJCut(double C6, double C12, double rcut) except +
    cdef cppclass  cLJCutPeriodic "pele::LJCutPeriodic":
        cLJCutPeriodic(double C6, double C12, double rcut, double * boxvec) except +
    cdef cppclass  cLJFrozen "pele::LJFrozen":
        cLJFrozen(double C6, double C12, _pele.Array[double] & reference_coords,
                  _pele.Array[np.int64_t] & frozen_dof) except +
    cdef cppclass  cLJ_Ilist "pele::LJ_interaction_list":
        cLJ_Ilist(_pele.Array[np.int64_t] & ilist, double C6, double C12) except +

cdef class LJ(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, eps=1.0, sigma=1.0, boxvec=None):
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            self.periodic = False
            self.thisptr = <_pele.cBasePotential*>new cLJ(4.*eps*sigma**6, 4.*eps*sigma**12)
        else:
            self.periodic = True
            bv = np.array(boxvec)
            self.thisptr = <_pele.cBasePotential*>new cLJPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12,
                                                                  <double*> bv.data)

cdef class LJCut(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, eps=1.0, sigma=1.0, rcut=2.5, boxvec=None):
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            self.periodic = False
            self.thisptr = <_pele.cBasePotential*>new cLJCut(4.*eps*sigma**6, 4.*eps*sigma**12, rcut)
        else:
            self.periodic = True
            bv = np.array(boxvec)
            self.thisptr = <_pele.cBasePotential*>new cLJCutPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12, rcut,
                                                                     <double*> bv.data)


cdef class LJFrozen(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, np.ndarray[double, ndim=1] reference_coords, 
                   frozen_atoms, 
                   eps=1.0, sigma=1.0, boxvec=None):
#        cdef np.ndarray[double, ndim=1] bv
        cdef np.ndarray[np.int64_t, ndim=1] frozen_dof
        frozen_dof = np.array([range(3*i,3*i+3) for i in frozen_atoms], dtype=int).reshape(-1)

        if boxvec is None:
            self.periodic = False
            self.thisptr = <_pele.cBasePotential*>new cLJFrozen(4.*eps*sigma**6, 4.*eps*sigma**12,
                        _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                        _pele.Array[np.int64_t](<np.int64_t *> frozen_dof.data, frozen_dof.size)
                        )
        else:
            self.periodic = True
#            bv = np.array(boxvec)
#            self.thisptr = <_pele.cBasePotential*>new cLJPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12,
#                                                                  <double*> bv.data)


cdef class LJInteractionList(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    def __cinit__(self, np.ndarray[np.int64_t, ndim=1] ilist, eps=1.0, sigma=1.0):
        self.thisptr = <_pele.cBasePotential*>new cLJ_Ilist( _pele.Array[np.int64_t](<np.int64_t*> ilist.data, <int> ilist.size), 4.*eps*sigma**6, 4.*eps*sigma**12)

