#cython: boundscheck=False
#cython: wraparound=False

from libcpp cimport bool as cbool
cimport cython
cimport numpy as np
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np

cdef extern from "pele/lowest_eig_potential.h" namespace "pele":
    cdef cppclass cOrthogonalizeTranslational "pele::OrthogonalizeTranslational":
        cOrthogonalizeTranslational(size_t natoms, size_t bdim, double tol) except+
        void orthogonalize(_pele.Array[double]&, _pele.Array[double]&) except+
    cdef cppclass cLowestEigPotential "pele::LowestEigPotential":
        cLowestEigPotential(shared_ptr[_pele.cBasePotential], _pele.Array[double], size_t, double) except +
        void reset_coords(_pele.Array[double]) except+