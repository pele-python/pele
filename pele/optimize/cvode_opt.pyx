"""
# distutils: language = C++
# cython: language_level=3str
"""
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import sys

import numpy as np

from pele.potentials import _pele
from pele.potentials cimport _pele
from pele.optimize import Result
from pele.potentials._pythonpotential import as_cpp_potential


cimport numpy as np
cimport pele.optimize._pele_opt as _pele_opt
from pele.optimize._pele_opt cimport shared_ptr
cimport cython
from cpython cimport bool as cbool

from libcpp cimport bool
from libcpp.vector cimport vector


cdef extern from "pele/cvode.hpp" namespace "pele":
    cpdef enum HessianType:
        DENSE,
        ITERATIVE

cdef extern from "pele/cvode.hpp" namespace "pele":
    cdef cppclass cppCVODEBDFOptimizer "pele::CVODEBDFOptimizer":
        cppCVODEBDFOptimizer(shared_ptr[_pele.cBasePotential], _pele.Array[double],
                             double, double, double, HessianType, bool, bool, int, double, double) except +
        double get_nhev() except +
        vector[double] get_time_trajectory() except +
        vector[double] get_gradient_norm_trajectory() except +
        vector[double] get_distance_trajectory() except +
        vector[double] get_energy_trajectory() except +
        vector[double] get_costly_time_trajectory() except +
        vector[vector[double]] get_coordinate_trajectory() except +
        vector[vector[double]] get_gradient_trajectory() except +
        

cdef class _Cdef_CVODEBDFOptimizer_CPP(_pele_opt.GradientOptimizer):
    """This class is the python interface for the c++ CVODEBDFOptimizer implementation
       Not quite an optimizer in the usual sense. Solves for the trajectory taken by
       the steepest descent path using the sundials solver CVODE_BDF. the solver terminates
       when it gets close to a minimum, just like an optimizer
    """
    cdef _pele.BasePotential pot
    cdef HessianType hess_type
    cdef bool save_trajectory
    def __cinit__(self, potential, x0,  double tol=1e-4, double rtol = 1e-4, double atol = 1e-4, int nsteps=10000, HessianType hessian_type=HessianType.DENSE, bool use_newton_stop_criterion=False, bool save_trajectory=False, int iterations_before_save=1, double newton_tol = 1e-5, double offset_factor = 1e-1):
        potential = as_cpp_potential(potential, verbose=True)

        self.pot = potential
        self.save_trajectory = save_trajectory
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        hess_type = hessian_type
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                new cppCVODEBDFOptimizer(self.pot.thisptr,
                             _pele.Array[double](<double*> x0c.data, x0c.size),
                                tol, rtol, atol, hess_type, use_newton_stop_criterion, save_trajectory, iterations_before_save, newton_tol, offset_factor))
    
    def get_result(self):
        cdef cppCVODEBDFOptimizer* cvode_ptr = <cppCVODEBDFOptimizer*> self.thisptr.get()
        res = super(_Cdef_CVODEBDFOptimizer_CPP, self).get_result()
        res["nhev"] = float(cvode_ptr.get_nhev())
        if self.save_trajectory:
            res["time_trajectory"] = np.array(cvode_ptr.get_time_trajectory())
            res["gradient_norm_trajectory"] = np.array(cvode_ptr.get_gradient_norm_trajectory())
            res["distance_trajectory"] = np.array(cvode_ptr.get_distance_trajectory())
            res["energy_trajectory"] = np.array(cvode_ptr.get_energy_trajectory())
            # costly trajectories
            res["costly_time_trajectory"] = np.array(cvode_ptr.get_costly_time_trajectory())
            res["coordinate_trajectory"] = np.array(cvode_ptr.get_coordinate_trajectory())
            res["gradient_trajectory"] = np.array(cvode_ptr.get_gradient_trajectory())
            

            
        return res
    
    def __reduce__(self):
        return (_Cdef_CVODEBDFOptimizer_CPP, (self.pot, self.x0, self.tol, self.rtol, self.atol, self.nsteps, self.iterative, self.use_newton_stop_criterion))

class CVODEBDFOptimizer_CPP(_Cdef_CVODEBDFOptimizer_CPP):
    """Python interface for the C++ CVODEBDFOptimizer implementation.
    """
