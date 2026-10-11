# distutils: language = c++
# cython: language_level=3str
# distutils: sources = modified_fire.cpp
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import numpy as np

from pele.potentials import _pele
from pele.potentials cimport _pele
from pele.potentials._pythonpotential import as_cpp_potential

cimport numpy as np
cimport pele.optimize._pele_opt as _pele_opt
from pele.optimize._pele_opt cimport shared_ptr
cimport cython

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "pele/optimizer.hpp" namespace "pele":
    cpdef enum StopCriterionType:
        GRADIENT,
        STEPNORM,
        NEWTON

cdef extern from "pele/modified_fire.hpp" namespace "pele":
    cdef cppclass cppMODIFIED_FIRE "pele::MODIFIED_FIRE":
        cppMODIFIED_FIRE(shared_ptr[_pele.cBasePotential], _pele.Array[double], 
                         double, double, double, size_t, double, double, double, double, double, bool, bool, int, StopCriterionType) except +
        vector[double] get_time_trajectory() except +
        vector[double] get_gradient_norm_trajectory() except +
        vector[double] get_distance_trajectory() except +
        vector[double] get_energy_trajectory() except +
        vector[double] get_costly_time_trajectory() except +
        vector[vector[double]] get_coordinate_trajectory() except +
        vector[vector[double]] get_gradient_trajectory() except +

cdef class _Cdef_MODIFIED_FIRE_CPP(_pele_opt.GradientOptimizer):
    """This class is the python interface for the c++ MODIFIED_FIRE implementation
    """  
    cdef _pele.BasePotential pot
    cdef double[:] x0
    cdef float dtstart
    cdef float dtmax
    cdef float maxstep
    cdef int Nmin
    cdef double finc
    cdef double fdec
    cdef double fa
    cdef double astart
    cdef double tol
    cdef bool stepback
    cdef int iprint
    cdef int nsteps
    cdef int verbosity
    cdef bool save_trajectory
    cdef int iterations_before_save
    cdef StopCriterionType stop_criterion
    
    
    def __cinit__(self, x0, potential, double dtstart = 0.1, double dtmax = 1, double maxstep=0.5, size_t Nmin=5, double finc=1.1, 
                   double fdec=0.5, double fa=0.99, double astart=0.1, double tol=1e-3, bool stepback = True, 
                   int iprint=-1, energy=None, gradient=None, int nsteps=10000, int verbosity=0, events = None, save_trajectory=False, iterations_before_save=1, StopCriterionType stop_criterion_type=StopCriterionType.GRADIENT):
        
        self.dtstart = dtstart
        self.dtmax = dtmax
        self.maxstep = maxstep
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.fa = fa
        self.astart = astart
        self.tol = tol
        self.stepback = stepback
        self.iprint = iprint
        self.nsteps = nsteps
        self.verbosity = verbosity
        self.save_trajectory = save_trajectory
        self.iterations_before_save = iterations_before_save
        self.stop_criterion = stop_criterion_type
        
        potential = as_cpp_potential(potential, verbose=verbosity>0)
        
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.x0 = x0c
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                        new cppMODIFIED_FIRE(pot.thisptr, _pele.Array[double](<double*> x0c.data, x0c.size),
                                             dtstart, dtmax, maxstep, Nmin, finc, fdec, fa, astart, tol, stepback, save_trajectory, iterations_before_save, stop_criterion_type))
        
        self.thisptr.get().set_max_iter(nsteps)
        self.thisptr.get().set_verbosity(verbosity)
        self.thisptr.get().set_iprint(iprint)
        self.pot = pot
        cdef np.ndarray[double, ndim=1] g_  
        if energy is not None and gradient is not None:
            g_ = gradient
            self.thisptr.get().set_func_gradient(energy, _pele.Array[double](<double*> g_.data, g_.size))

        self.events = events
        if self.events is None: 
            self.events = []
        
    def get_result(self):
        cdef cppMODIFIED_FIRE* cppfire = <cppMODIFIED_FIRE*> self.thisptr.get()
        result = super(_Cdef_MODIFIED_FIRE_CPP, self).get_result()
        if self.save_trajectory:
            result["time_trajectory"] = np.array(cppfire.get_time_trajectory())
            result["gradient_norm_trajectory"] = np.array(cppfire.get_gradient_norm_trajectory())
            result["distance_trajectory"] = np.array(cppfire.get_distance_trajectory())
            result["energy_trajectory"] = np.array(cppfire.get_energy_trajectory())
            # costly trajectories
            result["costly_time_trajectory"] = np.array(cppfire.get_costly_time_trajectory())
            result["coordinate_trajectory"] = np.array(cppfire.get_coordinate_trajectory())
            result["gradient_trajectory"] = np.array(cppfire.get_gradient_trajectory())
        return result
            
    def __reduce__(self):
        # energy, gradient and events are set to None, change if necessary later
        return (_Cdef_MODIFIED_FIRE_CPP, (self.x0, self.pot, self.dtstart, self.dtmax, self.maxstep, self.Nmin, self.finc, self.fdec, self.fa, self.astart, self.tol, self.stepback, self.iprint, None, None, self.nsteps, self.verbosity, None))

class ModifiedFireCPP(_Cdef_MODIFIED_FIRE_CPP):
    """Python interface for the C++ MODIFIED_FIRE implementation.
    """
    
#     def reset(self, coords):
#         """do one iteration"""
#         _Cdef_MODIFIED_FIRE_CPP.reset(self, coords)
       
