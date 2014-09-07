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
from pele.potentials._pele cimport array_wrap_np, pele_array_to_np
from libcpp.vector cimport vector as stdvector

# use external c++ class
cdef extern from "pele/aatopology.h" namespace "pele":
    cdef cppclass cppRBTopology "pele::RBTopology":
        cppRBTopology() except +
        void add_site(_pele.Array[double] x,
                      _pele.Array[double] cog,
                      double M, double W,
                      _pele.Array[double] S,
                      ) except +
        void get_zero_modes(_pele.Array[double]  x,
            stdvector[_pele.Array[double] ] & zev) except +
    cdef cppclass  cppRBPotentialWrapper "pele::RBPotentialWrapper":
        cppRBPotentialWrapper(shared_ptr[_pele.cBasePotential], shared_ptr[cppRBTopology]) except +
#         void add_site(_pele.Array[double]) except +

cdef class _cdef_RBTopology(object):
    cdef shared_ptr[cppRBTopology] thisptr
    def __cinit__(self, python_topology):
        self.thisptr = shared_ptr[cppRBTopology](new cppRBTopology() )
        
        # add sites to the c++ topology
        cdef np.ndarray[double, ndim=1] apos
        for site in python_topology.sites:
            apos = np.asarray(site.atom_positions, dtype=float).reshape(-1)
            self.thisptr.get().add_site(array_wrap_np(apos),
                                        array_wrap_np(site.cog),
                                        float(site.M), float(site.W), 
                                        array_wrap_np(site.S.reshape(-1))
                                        )
        
    def get_zero_modes(self, x):
        cdef np.ndarray[double, ndim=1] xnp = np.asarray(x, dtype=float, order="C")
        cdef stdvector[_pele.Array[double] ] c_zevlist
        self.thisptr.get().get_zero_modes(array_wrap_np(xnp), c_zevlist)
        zevlist = []
        for i in xrange(6):
            zevlist.append(pele_array_to_np(c_zevlist[i]))
        return zevlist
        


cdef class _cdef_RBPotentialWrapper(BasePotential):
    """define the python interface to the c++ RBPotentialWrapper implementation
    """
    topology = None
    def __cinit__(self, python_topology, potential):
        # wrap the potential so it can be used in the c++ classes
        if not issubclass(potential.__class__, _pele.BasePotential):
            potential = _pythonpotential.CppPotentialWrapper(potential)
        cdef _pele.BasePotential pot = potential

        # set up the cpp topology
        self.topology = _cdef_RBTopology(python_topology)
        cdef _cdef_RBTopology c_top = self.topology

        cdef cppRBPotentialWrapper * rbpot = new cppRBPotentialWrapper(pot.thisptr, c_top.thisptr)
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