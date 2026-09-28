"""
# distutils: language = C++
"""
cimport numpy as np
from libc.stdint cimport uintptr_t

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr

# use external c++ class
cdef extern from "pele/lammps_pele.h" namespace "pele":
    cdef cppclass  cLAMMPSPotential "pele::LAMMPSPotential":
        cLAMMPSPotential(uintptr_t handle) except +


cdef class LAMMPSPotential(_pele.BasePotential):
    """define the python interface to the LAMMPS library
    """
    def __cinit__(self, lmp):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLAMMPSPotential(<uintptr_t>lmp.lmp.value) )
