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

from libcpp cimport bool, nullptr
from libcpp.cast cimport static_cast


cdef extern from "pele/optimizer.hpp" namespace "pele":
    cdef enum StopCriterionType:
        GRADIENT,
        STEPNORM,
        NEWTON


cdef extern from "pele/extended_mixed_descent.hpp" namespace "pele":
    cdef cppclass cppExtendedMixedOptimizer "pele::ExtendedMixedOptimizer":
        cppExtendedMixedOptimizer(shared_ptr[_pele.cBasePotential], _pele.Array[double],
                                  shared_ptr[_pele.cBasePotential],
                                  double, int,
                                  double, double,
                                  double, double, bool, _pele.Array[double], bool, int, StopCriterionType) except +
        int get_nhev() except +
        int get_n_phase_1_steps() except +
        int get_n_phase_2_steps() except +
        int get_n_failed_phase_2_steps() except +

cdef class _Cdef_ExtendedMixedOptimizer_CPP(_pele_opt.GradientOptimizer):
    """Python interface for the C++ ExtendedMixedOptimizer implementation.
    """
    cdef _pele.BasePotential pot
    cdef _pele.BasePotential pot_ext
    cdef shared_ptr[_pele.cBasePotential] pot_ext_ptr
    cdef StopCriterionType stop_criterion_type
    cdef bool save_trajectory
    def __cinit__(self, potential, x0, potential_extension=None, global_symmetry_offset=[], double tol=1e-5,
                  int T=10, double step=1, double conv_tol = 1e-1,
                  double rtol=1e-3, double atol=1e-3, int nsteps=10000, bool iterative=False, bool save_trajectory=False, int steps_before_save=1, StopCriterionType stop_criterion_type=StopCriterionType.GRADIENT):
        potential = as_cpp_potential(potential, verbose=True)
        global_symmetry_offset = np.array(global_symmetry_offset, dtype=float)
        if len(global_symmetry_offset.shape) == 2:
            global_symmetry_offset = global_symmetry_offset.flatten()
        cdef np.ndarray[double, ndim=1] global_symmetry_offset_c = np.array(global_symmetry_offset, dtype=float)
        self.pot = potential
        self.stop_criterion_type = stop_criterion_type
        self.save_trajectory = save_trajectory
        if potential_extension is not None:
            pot_ext = as_cpp_potential(potential_extension, verbose=True)
            self.pot_ext = pot_ext
            self.pot_ext_ptr = self.pot_ext.thisptr
        else:
            self.pot_ext_ptr = static_cast[shared_ptr[_pele.cBasePotential]](nullptr)
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                new cppExtendedMixedOptimizer(self.pot.thisptr,
                                              _pele.Array[double](<double*> x0c.data, x0c.size), 
                                              self.pot_ext_ptr,
                                              tol, T, step, conv_tol, rtol, atol, iterative,
                                              _pele.Array[double](<double*> global_symmetry_offset_c.data,
                                                                  global_symmetry_offset_c.size), self.save_trajectory, steps_before_save, self.stop_criterion_type))
        cdef cppExtendedMixedOptimizer* mxopt_ptr = <cppExtendedMixedOptimizer*> self.thisptr.get()
    
    def get_result(self):
        cdef cppExtendedMixedOptimizer* mxopt_ptr = <cppExtendedMixedOptimizer*> self.thisptr.get()
        res = super(_Cdef_ExtendedMixedOptimizer_CPP, self).get_result()
        res["nhev"] = int(mxopt_ptr.get_nhev())
        res["n_phase_1"] = int(mxopt_ptr.get_n_phase_1_steps())
        res["n_phase_2"] = int(mxopt_ptr.get_n_phase_2_steps())
        res["n_failed_phase_2"] = int(mxopt_ptr.get_n_failed_phase_2_steps())
        return res

class ExtendedMixedOptimizer_CPP(_Cdef_ExtendedMixedOptimizer_CPP):
    """Python interface for the C++ ExtendedMixedOptimizer implementation.
    """
