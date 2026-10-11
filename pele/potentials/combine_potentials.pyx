"""
# distutils: language = C++
# cython: language_level=3str
"""

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
## combines potentials

cdef extern from "pele/combine_potentials.hpp" namespace "pele":
    cdef cppclass cppCombinedPotential "pele::CombinedPotential":
        cppCombinedPotential() except +
        void add_potential(shared_ptr[_pele.cBasePotential]) except +
        
        
        
        
        
        
        
cdef class CombinedPotential(_pele.BasePotential):
    """Sum native potentials added with ``add_potential``.

    Create an empty potential, then add each component before evaluating it.
    """
    def __cinit__(self):
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*>new cppCombinedPotential())
    def add_potential(self, _pele.BasePotential potential):
        (<cppCombinedPotential*>self.thisptr.get()).add_potential(potential.thisptr)
