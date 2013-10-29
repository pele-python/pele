cimport pele.potentials._pele as _pele
import numpy as np
cimport numpy as np
from cpython cimport bool

# use external c++ class
cdef extern from "lj.h" namespace "pele":
    cdef cppclass  cLJ "pele::LJ":
        cLJ(double C6, double C12) except +
    cdef cppclass  cLJPeriodic "pele::LJPeriodic":
        cLJPeriodic(double C6, double C12, double * boxvec) except +
    cdef cppclass  cLJFrozen "pele::LJFrozen":
        cLJFrozen(double C6, double C12, _pele.Array[double] & reference_coords,
                  _pele.Array[long] & frozen_dof) except +
    cdef cppclass  cLJ_Ilist "pele::LJ_interaction_list":
        cLJ_Ilist(_pele.Array[long] & ilist, double C6, double C12) except +
#    cdef cppclass  cLJAtomlist "pele::LJAtomList":
#        cLJAtomlist(double C6, double C12, 
#                    _pele.Array[size_t] & atoms1,
#                    _pele.Array[size_t] & atoms2,
#                    ) except +

cdef extern from "lj_cut.h" namespace "pele":
    cdef cppclass  cLJCut "pele::LJCut":
        cLJCut(double C6, double C12, double rcut) except +
    cdef cppclass  cLJCutPeriodic "pele::LJCutPeriodic":
        cLJCutPeriodic(double C6, double C12, double rcut, double * boxvec) except +
    cdef cppclass  cLJCutAtomlist "pele::LJCutAtomList":
        cLJCutAtomlist(double C6, double C12, double rcut,
                    _pele.Array[size_t] & atoms1,
                    _pele.Array[size_t] & atoms2) except +
        cLJCutAtomlist(double C6, double C12, double rcut,
                    _pele.Array[size_t] & atoms1) except +

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
        cdef np.ndarray[long, ndim=1] frozen_dof
        frozen_dof = np.array([range(3*i,3*i+3) for i in frozen_atoms], dtype=int).reshape(-1)

        if boxvec is None:
            self.periodic = False
            self.thisptr = <_pele.cBasePotential*>new cLJFrozen(4.*eps*sigma**6, 4.*eps*sigma**12,
                        _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                        _pele.Array[long](<long *> frozen_dof.data, frozen_dof.size)
                        )
        else:
            self.periodic = True
#            bv = np.array(boxvec)
#            self.thisptr = <_pele.cBasePotential*>new cLJPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12,
#                                                                  <double*> bv.data)


cdef class LJInteractionList(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    def __cinit__(self, np.ndarray[long, ndim=1] ilist, eps=1.0, sigma=1.0):
        self.thisptr = <_pele.cBasePotential*>new cLJ_Ilist( _pele.Array[long](<long*> ilist.data, <int> ilist.size), 4.*eps*sigma**6, 4.*eps*sigma**12)

cdef class BLJCut(_pele.BasePotential):
    def __cinit__(self, natoms, ntypeA, boxl=None, rcut=2.5, epsBB=0.5, sigBB=0.88, epsAB=1.5, sigAB=0.8):
        sigAA = 1.
        epsAA = 1.
        cdef np.ndarray[size_t, ndim=1] atomsAnp = np.array([range(ntypeA)],         dtype=np.uint).reshape(-1)
        cdef np.ndarray[size_t, ndim=1] atomsBnp = np.array([range(ntypeA, natoms)], dtype=np.uint).reshape(-1)
        cdef _pele.Array[size_t] atomsA = _pele.Array[size_t](<size_t*> atomsAnp.data, <size_t>atomsAnp.size)
        cdef _pele.Array[size_t] atomsB = _pele.Array[size_t](<size_t*> atomsBnp.data, <size_t>atomsBnp.size)
        cdef _pele.cCombinedPotential* combpot = new _pele.cCombinedPotential()
        combpot.add_potential(<_pele.cBasePotential*>new 
                                   cLJCutAtomlist(4.*epsAA*sigAA**6, 4.*epsAA*sigAA**12, rcut*sigAA,
                                                  atomsA
                                                  ))
        combpot.add_potential(<_pele.cBasePotential*>new 
                                   cLJCutAtomlist(4.*epsAB*sigAB**6, 4.*epsAB*sigAB**12, rcut*sigAB,
                                                  atomsA, atomsB
                                                  ))
        combpot.add_potential(<_pele.cBasePotential*>new 
                                   cLJCutAtomlist(4.*epsBB*sigBB**6, 4.*epsBB*sigBB**12, rcut*sigBB,
                                                  atomsB
                                                  ))
        self.thisptr = <_pele.cBasePotential*> combpot

    
    
