"""
# distutils: language = C++

This provides an interface to the c++ angle axis + rigid body classes.
Not everything is implemented in c++, just the parts that were bottlenecks
in python.
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
from libcpp cimport bool as cbool

# use external c++ class
cdef extern from "pele/aatopology.h" namespace "pele":

    cdef cppclass cppRigidFragment "pele::RigidFragment":
        cppRigidFragment(_pele.Array[double] x,
                      _pele.Array[double] cog,
                      double M, double W,
                      _pele.Array[double] S,
                      _pele.Array[double] inversion,
                      cbool can_invert
                      ) except +
        void add_symmetry_rotation(_pele.Array[double]) except +

    cdef cppclass cppRBTopology "pele::RBTopology":
        cppRBTopology() except +
        void add_site(cppRigidFragment & site) except +
        void get_zero_modes(_pele.Array[double]  x,
            stdvector[_pele.Array[double] ] & zev) except +
        double distance_squared(_pele.Array[double]  x1,
                              _pele.Array[double]  x2) except +
        # sn402: Overloaded
        double distance_squared_bulk(_pele.Array[double]  x1,
                              _pele.Array[double]  x2, 
                              _pele.Array[double] boxvec) except +
        void distance_squared_grad(_pele.Array[double]  x1,
                              _pele.Array[double]  x2,
                              _pele.Array[double]  grad) except +  
        # sn402: Overloaded                                                                                          
        void distance_squared_grad_bulk(_pele.Array[double]  x1,
                              _pele.Array[double]  x2,                            
                              _pele.Array[double]  grad,
                              _pele.Array[double] boxvec) except +    
                              
    cdef cppclass cppTransformAACluster "pele::TransformAACluster":
        cppTransformAACluster(cppRBTopology * top) except +
        void rotate(_pele.Array[double]  x,
                    _pele.Array[double]  mx) except +

    cdef cppclass cppMeasureAngleAxisCluster "pele::MeasureAngleAxisCluster":
        cppMeasureAngleAxisCluster(cppRBTopology * top) except +
        void align(_pele.Array[double]  x1,
                   _pele.Array[double]  x2) except +

    cdef cppclass  cppRBPotentialWrapper "pele::RBPotentialWrapper":
        cppRBPotentialWrapper(shared_ptr[_pele.cBasePotential], shared_ptr[cppRBTopology]) except +
#         void add_site(_pele.Array[double]) except +

cdef class _cppBaseTopology(object):
    pass

cdef class _cdef_RBTopology(_cppBaseTopology):
    cdef shared_ptr[cppRBTopology] thisptr
    def __cinit__(self, python_topology):
        self.thisptr = shared_ptr[cppRBTopology](new cppRBTopology() )
        
        for site in python_topology.sites:
            self._add_site(site)
    
    def _add_site(self, site):
        """ """
        cdef np.ndarray[double, ndim=1] apos = np.asarray(site.atom_positions, dtype=float).reshape(-1)
        S = np.asarray(site.S, dtype=float).reshape(-1)
        cdef cbool can_invert = site.inversion is not None
        if can_invert:
            inversion = np.asarray(site.inversion, dtype=float).reshape(-1)
            assert inversion.size == 9
        else:
            inversion = np.eye(3,3)
        cdef cppRigidFragment * rf = new cppRigidFragment(array_wrap_np(apos),
                                                          array_wrap_np(site.cog),
                                                          float(site.M), float(site.W), 
                                                          array_wrap_np(S),
                                                          array_wrap_np(inversion),
                                                          can_invert)
        for mx in site.symmetries:
            rf.add_symmetry_rotation(array_wrap_np(mx.reshape(-1)))
        self.thisptr.get().add_site(rf[0])
        del rf;
        
    def get_zero_modes(self, x):
        cdef np.ndarray[double, ndim=1] xnp = np.asarray(x, dtype=float, order="C")
        cdef stdvector[_pele.Array[double] ] c_zevlist
        self.thisptr.get().get_zero_modes(array_wrap_np(xnp), c_zevlist)
        zevlist = []
        for i in xrange(6):
            zevlist.append(pele_array_to_np(c_zevlist[i]))
        return zevlist
    
    def distance_squared(self, x1in, x2in):
        cdef np.ndarray[double, ndim=1] x1 = np.asarray(x1in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] x2 = np.asarray(x2in, dtype=float, order="C")
        cdef double d2 = self.thisptr.get().distance_squared(array_wrap_np(x1), array_wrap_np(x2))
        return d2
  
    # sn402: New    
    def distance_squared_bulk(self, x1in, x2in, boxvec):
        cdef np.ndarray[double, ndim=1] x1 = np.asarray(x1in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] x2 = np.asarray(x2in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] box = np.asarray(boxvec, dtype=float, order="C")
        cdef double d2 = self.thisptr.get().distance_squared_bulk(array_wrap_np(x1), array_wrap_np(x2), array_wrap_np(box))
        return d2
        
    def distance_squared_grad(self, x1in, x2in):
        cdef np.ndarray[double, ndim=1] x1 = np.asarray(x1in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] x2 = np.asarray(x2in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x1.size)
        self.thisptr.get().distance_squared_grad(array_wrap_np(x1), array_wrap_np(x2), array_wrap_np(grad))
        return grad

    # sn402: New    
    def distance_squared_grad_bulk(self, x1in, x2in, boxvec):
        cdef np.ndarray[double, ndim=1] x1 = np.asarray(x1in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] x2 = np.asarray(x2in, dtype=float, order="C")
        cdef np.ndarray[double, ndim=1] box = np.asarray(boxvec, dtype=float, order="C")           
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x1.size)

        self.thisptr.get().distance_squared_grad_bulk(array_wrap_np(x1), array_wrap_np(x2), array_wrap_np(grad), array_wrap_np(box))
        return grad    
        
cdef class _cdef_MeasureAngleAxisCluster(object):
    cdef shared_ptr[cppMeasureAngleAxisCluster] thisptr
    cdef _cdef_RBTopology topology
    def __cinit__(self, topology):
        # set up the cpp topology
        try:
            self.topology = topology
        except TypeError, e1:
            try:
                self.topology = topology.cpp_topology
            except (TypeError, AttributeError):
                print "can't set up the c++ topology"
                print "first tried", e1
                raise

        self.thisptr = shared_ptr[cppMeasureAngleAxisCluster](
                  new cppMeasureAngleAxisCluster(self.topology.thisptr.get()))
    
    def align(self, x1, x2):
        self.thisptr.get().align(array_wrap_np(x1), array_wrap_np(x1))

cdef class _cdef_TransformAACluster(object):
    cdef shared_ptr[cppTransformAACluster] thisptr
    cdef _cdef_RBTopology topology
    def __cinit__(self, topology):
        # set up the cpp topology
        try:
            self.topology = topology
        except TypeError, e1:
            try:
                self.topology = topology.cpp_topology
            except (TypeError, AttributeError):
                print "can't set up the c++ topology"
                print "first tried", e1
                raise

        self.thisptr = shared_ptr[cppTransformAACluster](
                  new cppTransformAACluster(self.topology.thisptr.get()))
    
    def rotate(self, x, mx):
        assert mx.size == 9
        self.thisptr.get().rotate(array_wrap_np(x), array_wrap_np(mx.reshape(-1)))


cdef class _cdef_RBPotentialWrapper(BasePotential):
    """define the python interface to the c++ RBPotentialWrapper implementation
    """
    cdef _cdef_RBTopology topology
    def __cinit__(self, topology, potential):
        # wrap the potential so it can be used in the c++ classes
        if not issubclass(potential.__class__, _pele.BasePotential):
            potential = _pythonpotential.CppPotentialWrapper(potential)
        cdef _pele.BasePotential pot = potential

        # set up the cpp topology
        try:
            self.topology = topology
        except TypeError, e1:
            try:
                self.topology = topology.cpp_topology
            except (TypeError, AttributeError):
                print "can't set up the c++ topology"
                print "first tried", e1
                raise
        
        assert self.topology is not None
            
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

class cdefTransformAACluster(_cdef_TransformAACluster):
    """Routines that apply transformations to a cluster of rigid bodies
    
    Parameters
    ----------
    topology : object
    """

class cdefMeasureAngleAxisCluster(_cdef_MeasureAngleAxisCluster):
    """Routines that perform measurements on a cluster of rigid bodies
    
    Parameters
    ----------
    topology : object
    """

class cdefRBTopology(_cdef_RBTopology):
    """Routines that perform measurements on a cluster of rigid bodies
    
    Parameters
    ----------
    python_topology : object
        This is the pythonic topology class which contains all the information 
        needed to construct the c++ topology class.
    """

