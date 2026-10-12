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

from libcpp cimport nullptr




cdef extern from "pele/generic_mixed_descent.hpp" namespace "pele":
    cdef cppclass cppGenericMixedDescent "pele::GenericMixedDescent":
        cppGenericMixedDescent(shared_ptr[_pele.cBasePotential], _pele.Array[double], 
                               double tol,
                               shared_ptr[_pele_opt.cGradientOptimizer],
                               shared_ptr[_pele_opt.cGradientOptimizer],
                               double, int, shared_ptr[_pele.cBasePotential]) except +
        cppGenericMixedDescent(shared_ptr[_pele.cBasePotential], _pele.Array[double], 
                               double tol,
                               shared_ptr[_pele_opt.cGradientOptimizer],
                               shared_ptr[_pele_opt.cGradientOptimizer],
                               double, int) except +
        double get_nhev() except +
        int get_n_phase_1_steps() except +
        int get_n_phase_2_steps() except +




cdef class _Cdef_GenericMixedDescent_CPP(_pele_opt.GradientOptimizer):
    """This class is the python interface for the c++ GenericMixedDescent implementation
    """
    cdef _pele.BasePotential pot
    cdef _pele.BasePotential pot_ext
    cdef _pele_opt.GradientOptimizer opt_conv
    cdef _pele_opt.GradientOptimizer opt_non_conv
    def __cinit__(self, potential, x0, optimizer_non_convex,
                  optimizer_convex, double tol, double translation_offset,
                  int steps_before_convex_check=1, potential_extension=None):
        potential = as_cpp_potential(potential, verbose=True)
        self.pot = potential
        if potential_extension is not None:
            self.pot_ext = as_cpp_potential(potential_extension, verbose=True)
        self.opt_conv = optimizer_convex
        self.opt_non_conv = optimizer_non_convex
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        
        # deals with python extension problems
        if potential_extension is not None:
            self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                new cppGenericMixedDescent(self.pot.thisptr,
                             _pele.Array[double](<double*> x0c.data, x0c.size), 
                             tol,
                             self.opt_non_conv.thisptr,
                             self.opt_conv.thisptr,
                             translation_offset,
                             steps_before_convex_check,
                             self.pot_ext.thisptr
                             ))
        else:
            self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                new cppGenericMixedDescent(self.pot.thisptr,
                             _pele.Array[double](<double*> x0c.data, x0c.size), 
                             tol,
                             self.opt_non_conv.thisptr,
                             self.opt_conv.thisptr,
                             translation_offset,
                             steps_before_convex_check,
                             ))
        cdef cppGenericMixedDescent* mxopt_ptr = <cppGenericMixedDescent*> self.thisptr.get()

    def get_result(self):
        cdef cppGenericMixedDescent* mxopt_ptr = <cppGenericMixedDescent*> self.thisptr.get()
        res = super(_Cdef_GenericMixedDescent_CPP, self).get_result()
        res["nhev"] = float(mxopt_ptr.get_nhev())
        res["n_phase_1"] = int(mxopt_ptr.get_n_phase_1_steps())
        res["n_phase_2"] = int(mxopt_ptr.get_n_phase_2_steps())
        return res

class GenericMixedDescent_CPP(_Cdef_GenericMixedDescent_CPP):
    """Python interface for the C++ GenericMixedDescent implementation.
    """
