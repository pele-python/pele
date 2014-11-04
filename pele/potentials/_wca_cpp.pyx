"""
# distutils: language = C++
"""
import numpy as np
from ctypes import c_size_t as size_t

cimport numpy as np
from cpython cimport bool

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr, array_wrap_np

# use external c++ class
cdef extern from "pele/wca.h" namespace "pele":
    cdef cppclass  cWCA "pele::WCA":
        cWCA(double sig, double eps) except +
    cdef cppclass  cWCAPeriodic "pele::WCAPeriodic":
        cWCAPeriodic(double sig, double eps, _pele.Array[double] boxvec) except +
    cdef cppclass  cWCANeighborList "pele::WCANeighborList":
        cWCANeighborList(_pele.Array[size_t] & ilist, double sig, double eps) except +
    cdef cppclass  cWCAAtomList "pele::WCAAtomList":
        cWCAAtomList(double sig, double eps, _pele.Array[size_t] atomlist) except +
    cdef cppclass  cWCA2D "pele::WCA2D":
        cWCA2D(double sig, double eps) except +
    cdef cppclass  cWCAPeriodic2D "pele::WCAPeriodic2D":
        cWCAPeriodic2D(double sig, double eps, _pele.Array[double] boxvec) except +

cdef class WCA(_pele.BasePotential):
    """define the python interface to the c++ WCA implementation
    """
    cpdef bool periodic 
    def __cinit__(self, sig=1.0, eps=1.0, ndim=3, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * ndim
        cdef np.ndarray[double, ndim=1] bv
        if boxvec is None:
            self.periodic = False
            if ndim == 2:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cWCA2D(sig, eps) )
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cWCA(sig, eps) )
        else:
            self.periodic = True
            assert(len(boxvec)==ndim)
            bv = np.array(boxvec, dtype=float)
            if ndim == 2:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                         cWCAPeriodic2D(sig, eps, _pele.Array[double](<double*> bv.data, bv.size)) )
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                         cWCAPeriodic(sig, eps, _pele.Array[double](<double*> bv.data, bv.size)) )
            
cdef class WCANeighborList(_pele.BasePotential):
    """define the python interface to the c++ WCA implementation
    """
    def __cinit__(self, ilist, eps=1.0, sigma=1.0):
        cdef np.ndarray[size_t, ndim=1] np_ilist = np.asarray(ilist, dtype=size_t).ravel()
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>
                        new cWCANeighborList( _pele.Array[size_t](<size_t*> np_ilist.data, <size_t> np_ilist.size),
                                                                     sigma, eps) )

cdef class WCAAtomList(_pele.BasePotential):
    """define the python interface to the c++ WCA implementation
    """
    def __cinit__(self, atoms, eps=1.0, sigma=1.0):
        cdef np.ndarray[size_t, ndim=1] atoms1  = np.array(atoms.reshape(-1), dtype=size_t) 

        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>
                        new cWCAAtomList(sigma, eps, 
                                         _pele.Array[size_t](<size_t*> atoms1.data, <int> atoms1.size)))

