"""
# distutils: language = C++
# cython: language_level=3str
"""
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import numpy as np
from ctypes import c_size_t as size_t

from pele.potentials import FrozenPotentialWrapper

cimport numpy as np
from cpython cimport bool

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np, array_wrap_np_size_t
# note: it is required to explicitly import BasePotential.  The compilation
# will fail if you try to use it as _pele.BasePotential.  I don't know why this
# is

# use external c++ class
cdef extern from "pele/frenkel.hpp" namespace "pele":
    cdef cppclass  cFrenkel "pele::Frenkel":
        cFrenkel(double sig, double eps, double rcut) except +
    cdef cppclass  cFrenkelPeriodic "pele::FrenkelPeriodic":
        cFrenkelPeriodic(double sig, double eps, double rcut, _pele.Array[double] boxvec) except +
    cdef cppclass cFrenkelPeriodicCellLists "pele::FrenkelPeriodicCellLists<3>":
        cFrenkelPeriodicCellLists(double sig, double eps, double rcut, _pele.Array[double] boxvec, double ncellx_scale)




cdef class Frenkel(_pele.BasePotential):
    """Python interface to the native Frenkel pair potential.
    """
    def __cinit__(self, natoms, sig=1.0, eps=1.0, rcut=2.0, boxvec=None, celllists=False, ncellx_scale=1.0):
        cdef np.ndarray[double, ndim=1] bv
        if boxvec==None and celllists==False:
            self.thisptr = _pele.shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cFrenkel(sig, eps, rcut))
        elif celllists==False:
            bv = np.array(boxvec)
            self.thisptr = _pele.shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cFrenkelPeriodic(sig, eps, rcut, array_wrap_np(bv)))
        elif celllists==True:
            bv = np.array(boxvec)
            self.thisptr = _pele.shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cFrenkelPeriodicCellLists(sig, eps, rcut, array_wrap_np(bv), ncellx_scale))



# cdef class Frenkelcpp(_pele.BasePotential):
#     """define the python interface to the c++ Frenkel implementation
#     """
#     def __cinit__(self, natoms, sig=1.0, eps=1.0, rcut=2.0):
#         self.thisptr = _pele.shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cFrenkel(sig, eps, rcut))
