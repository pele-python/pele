"""
# distutils: language = C++
"""
from ctypes import c_size_t as size_t

cimport numpy as np
import numpy as np

from pele.potentials._pele cimport shared_ptr, array_wrap_np
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport BasePotential


# define the external c++ classes
cdef extern from "pele/wca.h" namespace "pele":
    cdef cppclass  cWCAAtomList "pele::WCAAtomList":
        cWCAAtomList(double sig, double eps, _pele.Array[size_t] atomlist) except +
    cdef cppclass  cWCA "pele::WCA":
        cWCA(double sig, double eps) except +


cdef extern from "pele/lj_cut.h" namespace "pele":
    cdef cppclass  cLJCutAtomlist "pele::LJCutAtomList":
        cLJCutAtomlist(double C6, double C12, double rcut,
                    _pele.Array[size_t] & atoms1,
                    _pele.Array[size_t] & atoms2) except +
        cLJCutAtomlist(double C6, double C12, double rcut,
                    _pele.Array[size_t] & atoms1) except +

cdef extern from "pele/harmonic.h" namespace "pele":
#    cdef cppclass  cHarmonicAtomList "pele::HarmonicAtomList":
#        cHarmonicAtomList(double k,
#                    _pele.Array[size_t] & atoms1,
#                    _pele.Array[size_t] & atoms2) except +
#        cHarmonicAtomList(double k,
#                    _pele.Array[size_t] & atoms1) except +
    cdef cppclass  cHarmonicNeighborList "pele::HarmonicNeighborList":
        cHarmonicNeighborList(double k, _pele.Array[size_t] & ilist) except +



cdef class PlatePotential(BasePotential):
    """define the python interface to the c++ implementation
    """
    def __cinit__(self, harmonic_atoms1, harmonic_atoms2, lj_atoms, k=1.):
        assert harmonic_atoms1.size == harmonic_atoms2.size
        harmonic_nlist = np.zeros([harmonic_atoms1.size, 2])
        harmonic_nlist[:,0] = harmonic_atoms1
        harmonic_nlist[:,1] = harmonic_atoms2
        
        
        cdef np.ndarray[size_t, ndim=1] np_harmonic_nlist = np.array(harmonic_nlist.reshape(-1), dtype=size_t) 
#        cdef np.ndarray[size_t, ndim=1] np_harmonic_atoms2 = np.array(harmonic_atoms2.reshape(-1), dtype=size_t) 
        cdef np.ndarray[size_t, ndim=1] np_lj_atoms = np.array(lj_atoms.reshape(-1), dtype=size_t) 
#        cdef np.ndarray[size_t, ndim=1] np_wca_atoms = np.array(wca_atoms.reshape(-1), dtype=size_t) 

        # Wrap them with pele arrays
        cdef _pele.Array[size_t] p_harmonic_nlist = _pele.Array[size_t](<size_t*> np_harmonic_nlist.data, 
                                                                         <size_t>  np_harmonic_nlist.size)
#        cdef _pele.Array[size_t] p_harmonic_atoms2 = _pele.Array[size_t](<size_t*> np_harmonic_atoms2.data, 
#                                                                         <size_t>  np_harmonic_atoms2.size)
        cdef _pele.Array[size_t] p_lj_atoms = _pele.Array[size_t](<size_t*> np_lj_atoms.data, 
                                                                         <size_t>  np_lj_atoms.size)
#        cdef _pele.Array[size_t] p_wca_atoms = _pele.Array[size_t](<size_t*> np_wca_atoms.data, 
#                                                                         <size_t>  np_wca_atoms.size)


        # make the combined potential
        cdef _pele.cCombinedPotential* combpot = new _pele.cCombinedPotential()
        
        # add the potential for the harmonic
        combpot.add_potential(shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                   cHarmonicNeighborList(k, p_harmonic_nlist)
                                                  ))
        
        # add the potential for the lj
        cdef double c6 =  4.
        cdef double c12 = 4.
        cdef double rcut = 100.
        combpot.add_potential(shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                   cLJCutAtomlist(c6, c12, rcut, p_lj_atoms
                                                  )))
        
        # add the WCA potential
        combpot.add_potential(shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                                   cWCA(1., 1.)))

        
        # save the combined potential in the format of a shared_ptr as self.thisptr
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> combpot )

