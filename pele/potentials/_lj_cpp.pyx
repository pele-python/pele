"""
# distutils: language = C++
"""
import numpy as np
from ctypes import c_size_t as size_t

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
cimport numpy as np
from cpython cimport bool


# use external c++ class
cdef extern from "pele/lj.h" namespace "pele":
    cdef cppclass  cLJ "pele::LJ":
        cLJ(double C6, double C12) except +
    cdef cppclass  cLJPeriodic "pele::LJPeriodic":
        cLJPeriodic(double C6, double C12, _pele.Array[double] boxvec) except +
    cdef cppclass  cLJFrozen "pele::LJFrozen":
        cLJFrozen(double C6, double C12, _pele.Array[double] & reference_coords,
                  _pele.Array[size_t] & frozen_dof) except +
    cdef cppclass  cLJNeighborList "pele::LJNeighborList":
        cLJNeighborList(_pele.Array[long] & ilist, double C6, double C12) except +

cdef extern from "pele/lj_cut.h" namespace "pele":
    cdef cppclass  cLJCut "pele::LJCut":
        cLJCut(double C6, double C12, double rcut) except +
    cdef cppclass  cLJCutPeriodic "pele::LJCutPeriodic":
        cLJCutPeriodic(double C6, double C12, double rcut, _pele.Array[double] boxvec) except +
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
    def __cinit__(self, eps=1.0, sig=1.0, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * 3
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            self.periodic = False
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLJ(4.*eps*sig**6, 4.*eps*sig**12) )
        else:
            self.periodic = True
            bv = np.array(boxvec, dtype=float)
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLJPeriodic(4.*eps*sig**6, 4.*eps*sig**12,
                                                              _pele.Array[double](<double*> bv.data, bv.size)) )

cdef class LJCut(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, eps=1.0, sigma=1.0, rcut=2.5, boxvec=None):
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            self.periodic = False
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLJCut(4.*eps*sigma**6, 4.*eps*sigma**12, rcut) )
        else:
            self.periodic = True
            bv = np.array(boxvec)
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLJCutPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12, rcut,
                                                                 _pele.Array[double](<double*> bv.data, bv.size)) )


cdef class LJFrozen(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, np.ndarray[double, ndim=1] reference_coords, 
                   frozen_atoms, 
                   eps=1.0, sigma=1.0, boxvec=None):
#        cdef np.ndarray[double, ndim=1] bv
        cdef np.ndarray[size_t, ndim=1] frozen_dof
        frozen_dof = np.array([range(3*i,3*i+3) for i in frozen_atoms], dtype=int).reshape(-1)

        if boxvec is None:
            self.periodic = False
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLJFrozen(4.*eps*sigma**6, 4.*eps*sigma**12,
                        _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                        _pele.Array[size_t](<size_t *> frozen_dof.data, frozen_dof.size)
                        ) )
        else:
            self.periodic = True
            raise NotImplementedError("periodic LJ with frozen atoms is not implemented yet")
#            bv = np.array(boxvec)
#            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLJPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12,
#                                                                  <double*> bv.data) )


cdef class LJNeighborList(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    def __cinit__(self, np.ndarray[long, ndim=1] ilist, eps=1.0, sigma=1.0):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new cLJNeighborList( 
                     _pele.Array[long](<long*> ilist.data, <int> ilist.size), 4.*eps*sigma**6, 4.*eps*sigma**12) )

cdef class BLJCut(_pele.BasePotential):
    def __cinit__(self, natoms, ntypeA, boxl=None, rcut=2.5, epsAA=1., sigAA=1., 
                   epsBB=0.5, sigBB=0.88, epsAB=1.5, sigAB=0.8):
        cdef np.ndarray[size_t, ndim=1] atomsAnp = np.array([range(ntypeA)],         dtype=size_t).reshape(-1)
        cdef np.ndarray[size_t, ndim=1] atomsBnp = np.array([range(ntypeA, natoms)], dtype=size_t).reshape(-1)
        cdef _pele.Array[size_t] atomsA = _pele.Array[size_t](<size_t*> atomsAnp.data, <size_t>atomsAnp.size)
        cdef _pele.Array[size_t] atomsB = _pele.Array[size_t](<size_t*> atomsBnp.data, <size_t>atomsBnp.size)
        cdef _pele.cCombinedPotential* combpot = new _pele.cCombinedPotential()
        combpot.add_potential(shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                   cLJCutAtomlist(4.*epsAA*sigAA**6, 4.*epsAA*sigAA**12, rcut*sigAA,
                                                  atomsA
                                                  )))
        combpot.add_potential(shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                   cLJCutAtomlist(4.*epsAB*sigAB**6, 4.*epsAB*sigAB**12, rcut*sigAB,
                                                  atomsA, atomsB
                                                  )))
        combpot.add_potential(shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                   cLJCutAtomlist(4.*epsBB*sigBB**6, 4.*epsBB*sigBB**12, rcut*sigBB,
                                                  atomsB
                                                  )))
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> combpot )

    
cdef class _ErrorPotential(_pele.BasePotential):
    """this is a test potential which should raise an exception when called
    
    for testing purposes only
    """
    def __cinit__(self):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new _pele.cBasePotential() )
  
