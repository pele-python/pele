"""
# distutils: language = C++
"""
import numpy as np
from ctypes import c_size_t as size_t

from pele.potentials import FrozenPotentialWrapper

cimport numpy as np
from cpython cimport bool

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np, array_wrap_np_size_t


# use external c++ class
cdef extern from "pele/lj.h" namespace "pele":
    cdef cppclass  cLJ "pele::LJ":
        cLJ(double C6, double C12) except +
    cdef cppclass  cLJPeriodic "pele::LJPeriodic":
        cLJPeriodic(double C6, double C12, _pele.Array[double] boxvec) except +
    cdef cppclass  cLJNeighborList "pele::LJNeighborList":
        cLJNeighborList(_pele.Array[size_t] & ilist, double C6, double C12) except +

cdef extern from "pele/lj_cut.h" namespace "pele":
    cdef cppclass  cLJCut "pele::LJCut":
        cLJCut(double C6, double C12, double rcut) except +
    cdef cppclass  cLJCutPeriodic "pele::LJCutPeriodic":
        cLJCutPeriodic(double C6, double C12, double rcut, _pele.Array[double] boxvec) except +
    cdef cppclass  cLJCutAtomlist "pele::LJCutAtomList":
        cLJCutAtomlist(double C6, double C12, double rcut,
                    _pele.Array[size_t] atoms1,
                    _pele.Array[size_t] atoms2) except +
        cLJCutAtomlist(double C6, double C12, double rcut,
                    _pele.Array[size_t] atoms1) except +
    cdef cppclass  cLJCutPeriodicAtomList "pele::LJCutPeriodicAtomList":
        cLJCutPeriodicAtomList(double C6, double C12, double rcut, _pele.Array[double] boxvec,
                    _pele.Array[size_t] atoms1,
                    _pele.Array[size_t] atoms2) except +
        cLJCutPeriodicAtomList(double C6, double C12, double rcut, _pele.Array[double] boxvec,
                    _pele.Array[size_t] atoms1) except +
    cdef cppclass cppLJCutPeriodicCellLists "pele::LJCutPeriodicCellLists<3>":
        cppLJCutPeriodicCellLists(double C6, double C12, double rcut, 
                                  _pele.Array[double] boxvec, double ncellx_scale) except +

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
                                                              array_wrap_np(bv)) )

cdef class LJCut(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, eps=1.0, sigma=1.0, rcut=2.5, boxvec=None):
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            self.periodic = False
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                     cLJCut(4.*eps*sigma**6, 4.*eps*sigma**12, rcut) )
        else:
            self.periodic = True
            bv = np.array(boxvec)
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                     cLJCutPeriodic(4.*eps*sigma**6, 4.*eps*sigma**12, rcut,
                                    array_wrap_np(bv)) )

cdef class LJCutCellLists(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    cpdef bool periodic 
    def __cinit__(self, eps=1.0, sigma=1.0, rcut=2.5, boxvec=None, ncellx_scale=1.):
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            raise NotImplementedError("LJCutCellLists currently only works with periodic bounds")
            self.periodic = False
        else:
            self.periodic = True
            bv = np.array(boxvec)
            self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
                     cppLJCutPeriodicCellLists(4.*eps*sigma**6, 4.*eps*sigma**12, rcut,
                                               array_wrap_np(bv), ncellx_scale))


cdef class LJFrozen(_pele.BasePotential):
    """Lennard Jones potential with frozen atoms
    
    note: this should not really be used.  It is preferable to just create the potential
    and wrap it manually with FrozenPotentialWrapper
    """
    cpdef bool periodic 
    def __init__(self, np.ndarray[double, ndim=1] reference_coords, 
                   frozen_atoms, eps=1.0, sigma=1.0, boxvec=None):
        frozen_dof = np.array([range(3*i,3*i+3) for i in frozen_atoms], dtype=int).reshape(-1)

        lj = LJ(eps=eps, sigma=sigma, boxvec=boxvec)
        self.periodic = boxvec is not None
        cdef _pele.BasePotential frozen_pot_wrapper = FrozenPotentialWrapper(lj, reference_coords, frozen_dof)
        self.thisptr = frozen_pot_wrapper.thisptr
        

cdef class LJNeighborList(_pele.BasePotential):
    """define the python interface to the c++ LJ implementation
    """
    def __cinit__(self, ilist, eps=1.0, sigma=1.0):
        cdef np.ndarray[size_t, ndim=1] np_ilist = np.array(ilist, dtype=size_t).reshape(-1)
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new cLJNeighborList( 
                     _pele.Array[size_t](<size_t*> np_ilist.data, <size_t> np_ilist.size), 4.*eps*sigma**6, 4.*eps*sigma**12) )

cdef class BLJCut(_pele.BasePotential):
    """Binary Lennard-Jones with a cutoff
    """
    def __cinit__(self, natoms, ntypeA, boxvec=None, rcut=2.5, epsAA=1., sigAA=1., 
                   epsBB=0.5, sigBB=0.88, epsAB=1.5, sigAB=0.8):

#        # put the box lengths into a numpy array
#        cdef np.ndarray[double, ndim=1] bv = np.array([0., 0., 0.])
#        if boxvec is None:
#            periodic = False
#        else:
#            periodic = True
#            bv = np.array(boxvec, dtype=float)
        
        # make the lists of atom indices.  This could be passed
        cdef np.ndarray[size_t, ndim=1] atomsAnp = np.array(range(ntypeA),         dtype=size_t).reshape(-1)
        cdef np.ndarray[size_t, ndim=1] atomsBnp = np.array(range(ntypeA, natoms), dtype=size_t).reshape(-1)
        
        # make the combined potential
        cdef _pele.cCombinedPotential* combpot = new _pele.cCombinedPotential()
        
        # create the potentials
        ljAA = LJCutAtomList(atomsAnp, rcut=rcut, eps=epsAA, sig=sigAA, boxvec=boxvec)
        ljAB = LJCutAtomList(atomsAnp, atoms2=atomsBnp, rcut=rcut, eps=epsAB, sig=sigAB, boxvec=boxvec)
        ljBB = LJCutAtomList(atomsBnp, rcut=rcut, eps=epsBB, sig=sigBB, boxvec=boxvec)

        # add them to the combined pot
        combpot.add_potential(ljAA.thisptr)
        combpot.add_potential(ljAB.thisptr)
        combpot.add_potential(ljBB.thisptr)

        # save the combined potential in the format of a shared_ptr as self.thisptr
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> combpot )

cdef class LJCutAtomList(_pele.BasePotential):
    """define the python interface to the c++ implementation
    """
    def __cinit__(self, atoms1, atoms2=None, eps=1.0, sig=1.0, rcut=2.5, boxvec=None):
        cdef np.ndarray[size_t, ndim=1] atoms1_np  = np.array(atoms1.reshape(-1), dtype=size_t) 
        cdef _pele.Array[size_t] atoms1_a = _pele.Array[size_t](<size_t*> atoms1_np.data, <size_t>atoms1_np.size)

        cdef np.ndarray[size_t, ndim=1] atoms2_np 
        cdef _pele.Array[size_t] atoms2_a
        
        cdef np.ndarray[double, ndim=1] bv_np = np.array([0., 0., 0.])
        cdef _pele.Array[double] bv
        if boxvec is None:
            periodic = False
        else:
            periodic = True
            bv_np = np.array(boxvec, dtype=float)
            bv = array_wrap_np(bv_np)
            bv = bv.copy()

        
        if atoms2 is None:
            if periodic:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>
                                           new cLJCutPeriodicAtomList(4.*eps*sig**6, 4.*eps*sig**12, rcut*sig, bv,
                                                              atoms1_a))
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>
                                           new cLJCutAtomlist(4.*eps*sig**6, 4.*eps*sig**12, rcut*sig,
                                                              atoms1_a))
        else:
            atoms2_np  = np.array(atoms2.reshape(-1), dtype=size_t) 
            atoms2_a = _pele.Array[size_t](<size_t*> atoms2_np.data, <size_t>atoms2_np.size)
            if periodic:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>
                                           new cLJCutPeriodicAtomList(4.*eps*sig**6, 4.*eps*sig**12, rcut*sig, bv,
                                                              atoms1_a, atoms2_a))
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>
                                           new cLJCutAtomlist(4.*eps*sig**6, 4.*eps*sig**12, rcut*sig,
                                                              atoms1_a, atoms2_a))

    
cdef class _ErrorPotential(_pele.BasePotential):
    """this is a test potential which should raise an exception when called
    
    for testing purposes only
    """
    def __cinit__(self):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new _pele.cBasePotential() )
  
