"""
# distutils: language = C++

basic potential interface stuff    
"""

from ctypes import c_size_t as size_t

cimport numpy as np
import numpy as np

from pele.potentials._pythonpotential import as_cpp_potential
from pele.potentials cimport _pele
from pele.potentials._pele cimport Array, array_wrap_np, shared_ptr

#===============================================================================
# pele::FrozenPotentialWrapper
#===============================================================================
cdef extern from "pele/frozen_atoms.h" namespace "pele":
    cdef cppclass  cppFrozenPotentialWrapper "pele::FrozenPotentialWrapper":
        cppFrozenPotentialWrapper(shared_ptr[_pele.cBasePotential] potential,
            Array[double] reference_coords,
            Array[size_t] frozen_dof) except +
        Array[double] get_full_coords(Array[double] reduced_coords) except +
        Array[double] get_reduced_coords(Array[double] full_coords) except +
        size_t ndof() except +
        Array[size_t] get_mobile_dof() except +
        Array[size_t] get_frozen_dof() except +





cdef class FrozenPotentialWrapper(_pele.BasePotential):
    """Wrapper for a potential object for freezing degrees of freedom
    
    Parameters
    ----------
    potential : object
        the pele potential object to be wrapped
    reference_coords : numpy array
        a set of reference coordinates.  This defines the positions of the frozen 
        coordinates and the total number of degrees of freedom
    frozen_dof : list
        a list of the frozen degrees of freedom
        
    Notes
    -----
    
    This uses class FrozenCoordsConverter to convert back and forth between the full
    set of coordinates and the reduced set of coordinates (those that are still mobile).
    
    The functions getEnergy, etc, accept a reduced set of coordinates,
    and passes it to the wrapped potential for calculation of the energy.  
    
    You can convert between the full and reduced representation using the functions
    `get_reduced_coords()` and `get_full_coords()`.
    
    Examples
    --------
    The following example shows how to wrap the lennard jones potential and freeze
    the first 6 degrees of freedom (2 atoms).  It then does a minimization on the 
    reduced coordinates and prints off some information
    
        import numpy as np
        from pele.potentials import LJ, FrozenPotentialWrapper
        from pele.optimize import mylbfgs
        natoms = 4
        pot = LJ()
        
        reference_coords = np.random.uniform(-1, 1, [3*natoms])
        print reference_coords
        
        # freeze the first two atoms (6 degrees of freedom)
        frozen_dof = range(6)
        
        fpot = FrozenPotentialWrapper(pot, reference_coords, frozen_dof)
        
        reduced_coords = fpot.coords_converter.get_reduced_coords(reference_coords)
        
        print "the energy in the full representation:" 
        print pot.getEnergy(reference_coords)
        print "is the same as the energy in the reduced representation:"
        print fpot.getEnergy(reduced_coords)
        
        ret = mylbfgs(reduced_coords, fpot)
        print "after a minimization the energy is ", ret.energy, "and the rms gradient is", ret.rms
        print "the coordinates of the frozen degrees of freedom are unchanged"
        print "starting coords:", reference_coords
        print "minimized coords:", fpot.coords_converter.get_full_coords(ret.coords)

    """
    cdef cppFrozenPotentialWrapper *direct_ptr # direct pointer for convenience only
    def __init__(self, potential, reference_coords, frozen_dof):
        cdef _pele.BasePotential pot
        # wrap the potential if necessary 
        pot = as_cpp_potential(potential)
        
        cdef np.ndarray[double, ndim=1] reference_coords_np = np.asarray(reference_coords, dtype=float).ravel()
        cdef np.ndarray[size_t, ndim=1] frozen_dof_np = np.asarray(frozen_dof, dtype=size_t).ravel()
        
        
        self.direct_ptr = new cppFrozenPotentialWrapper(
                pot.thisptr, array_wrap_np(reference_coords_np), _pele.array_wrap_np_size_t(frozen_dof_np))
        # now the shared_ptr takes ownership of it
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> self.direct_ptr )

        
    def get_full_coords(self, reduced_coords):
        cdef Array[double] x = self.direct_ptr.get_full_coords(array_wrap_np(reduced_coords))
        return _pele.pele_array_to_np(x)
    
    def get_reduced_coords(self, full_coords):
        cdef Array[double] x = self.direct_ptr.get_reduced_coords(array_wrap_np(full_coords))
        return _pele.pele_array_to_np(x)
    
    def get_mobile_dof(self):
        return _pele.pele_array_to_np_size_t(self.direct_ptr.get_mobile_dof())
    
    def get_frozen_dof(self):
        return _pele.pele_array_to_np_size_t(self.direct_ptr.get_frozen_dof())



