"""
# distutils: language = C++
"""
import numpy as np
from pele.potentials import _pele
from pele.potentials import _pythonpotential

cimport numpy as np
from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport BasePotential



# use external c++ class
cdef extern from "pele/aatopology.h" namespace "pele":
    cdef cppclass cppRBTopology "pele::RBTopology":
        cppRBTopology() except +
        void add_site(_pele.Array[double]) except +
    cdef cppclass  cppRBPotentialWrapper "pele::RBPotentialWrapper":
        cppRBPotentialWrapper(shared_ptr[_pele.cBasePotential], shared_ptr[cppRBTopology]) except +
#         void add_site(_pele.Array[double]) except +

cdef class _cdef_RBTopology(object):
    cdef shared_ptr[cppRBTopology] thisptr
    def __cinit__(self, python_topology):
        self.thisptr = shared_ptr[cppRBTopology](new cppRBTopology() )
        cdef np.ndarray[double, ndim=1] apos
        for site in python_topology.sites:
            apos = np.asarray(site.atom_positions, dtype=float).reshape(-1)
            self.thisptr.get().add_site(_pele.Array[double] (<double *> apos.data, apos.size) )


cdef class _cdef_RBPotentialWrapper(BasePotential):
    """define the python interface to the c++ RBPotentialWrapper implementation
    """
    cdef _cdef_RBTopology topology
    def __cinit__(self, python_topology, potential):
        # wrap the potential so it can be used in the c++ classes
        if not issubclass(potential.__class__, _pele.BasePotential):
            potential = _pythonpotential.CppPotentialWrapper(potential)
        cdef _pele.BasePotential pot = potential

        # set up the cpp topology
        self.topology = _cdef_RBTopology(python_topology)

        cdef cppRBPotentialWrapper * rbpot = new cppRBPotentialWrapper(pot.thisptr, self.topology.thisptr)
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> rbpot )

class RBPotentialWrapper(_cdef_RBPotentialWrapper):
    """c++ potential wrapper for Rigid Body potentials
    
    Parameters
    -----------
    topology : class derived from RBTopology
        this defines all the information about the topology of the rigid body system,
        in particular how to convert from rigid body coords to atomistic coords and back.
    atomistic_potential : pele potential
        this class will compute the energy and gradient in atomistic coordinates
    """