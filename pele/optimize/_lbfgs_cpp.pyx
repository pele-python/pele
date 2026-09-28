"""
# distutils: language = C++
# cython: language_level=3str
"""
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import sys
from collections import namedtuple

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


cdef extern from "pele/optimizer.hpp" namespace "pele":
    cpdef enum StopCriterionType:
        GRADIENT,
        STEPNORM,
        NEWTON


# import the externally defined lbfgs implementation
cdef extern from "pele/lbfgs.hpp" namespace "pele":
    cdef cppclass cppLBFGS "pele::LBFGS":
        cppLBFGS(shared_ptr[_pele.cBasePotential], _pele.Array[double], double, int, bool, int, StopCriterionType) except +

        void set_H0(double) except +
        void set_tol(double) except +
        void set_maxstep(double) except +
        void set_max_f_rise(double) except +
        void set_use_relative_f(int) except +
        void set_max_iter(int) except +
        void set_iprint(int) except +
        void set_verbosity(int) except +

        double get_H0() except +
        _pele.Array[double] get_s() except +
        _pele.Array[double] get_y() except +
        _pele.Array[double] get_rho() except +
        int get_k() except +
        void set_state(_pele.Array[double], _pele.Array[double], _pele.Array[double], double, int) except +



LBFGSState = namedtuple("LBFGSState", "s y rho H0 k")


cdef class _Cdef_LBFGS_CPP(_pele_opt.GradientOptimizer):
    """This class is the python interface for the c++ LBFGS implementation
    """
    cdef _pele.BasePotential pot
    cdef StopCriterionType stop_criterion
    cdef bool save_trajectory
    
    def __cinit__(self, x0, potential, double tol=1e-5, int M=4, double maxstep=0.1, 
                  double maxErise=1e-4, double H0=0.1, int iprint=-1,
                  energy=None, gradient=None,
                  int nsteps=10000, int verbosity=0, events=None, logger=None,
                  rel_energy=False, bool save_trajectory=False, int iterations_before_save=1, StopCriterionType stop_criterion=StopCriterionType.GRADIENT):
        potential = as_cpp_potential(potential, verbose=verbosity>0)

        self.pot = potential
        self.save_trajectory = save_trajectory
        self.stop_criterion = stop_criterion
        if logger is not None:
            print("warning c++ LBFGS is ignoring logger")
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                new cppLBFGS(self.pot.thisptr, 
                             _pele.Array[double](<double*> x0c.data, x0c.size),
                             tol, M, self.save_trajectory, iterations_before_save, self.stop_criterion) )
        cdef cppLBFGS* lbfgs_ptr = <cppLBFGS*> self.thisptr.get()
        lbfgs_ptr.set_H0(H0)
        lbfgs_ptr.set_maxstep(maxstep)
        lbfgs_ptr.set_max_f_rise(maxErise)
        lbfgs_ptr.set_max_iter(nsteps)
        lbfgs_ptr.set_verbosity(verbosity)
        lbfgs_ptr.set_iprint(iprint)
        if rel_energy:
            lbfgs_ptr.set_use_relative_f(1)
        
        cdef np.ndarray[double, ndim=1] g_  
        if energy is not None and gradient is not None:
            g_ = gradient
            self.thisptr.get().set_func_gradient(energy, _pele.Array[double](<double*> g_.data, g_.size))

        self.events = events
        if self.events is None: 
            self.events = []
    
    def set_H0(self, H0):
        cdef cppLBFGS* lbfgs_ptr = <cppLBFGS*> self.thisptr.get()
        lbfgs_ptr.set_H0(float(H0))
    
    def get_state(self):
        """return the LBFGS memory, see set_state"""
        cdef cppLBFGS* lbfgs_ptr = <cppLBFGS*> self.thisptr.get()
        return LBFGSState(s=_pele.pele_array_to_np(lbfgs_ptr.get_s()),
                          y=_pele.pele_array_to_np(lbfgs_ptr.get_y()),
                          rho=_pele.pele_array_to_np(lbfgs_ptr.get_rho()),
                          H0=float(lbfgs_ptr.get_H0()),
                          k=int(lbfgs_ptr.get_k()))

    def set_state(self, state):
        """restore the LBFGS memory from get_state of an optimizer with the same size and M

        This keeps the curvature information when restarting at new coordinates.
        """
        cdef cppLBFGS* lbfgs_ptr = <cppLBFGS*> self.thisptr.get()
        cdef np.ndarray[double, ndim=1] s = np.ascontiguousarray(state.s, dtype=float).ravel()
        cdef np.ndarray[double, ndim=1] y = np.ascontiguousarray(state.y, dtype=float).ravel()
        cdef np.ndarray[double, ndim=1] rho = np.ascontiguousarray(state.rho, dtype=float).ravel()
        lbfgs_ptr.set_state(_pele.Array[double](<double*> s.data, s.size),
                            _pele.Array[double](<double*> y.data, y.size),
                            _pele.Array[double](<double*> rho.data, rho.size),
                            float(state.H0), int(state.k))

    def get_result(self):
        cdef cppLBFGS* lbfgs_ptr = <cppLBFGS*> self.thisptr.get()
        res = super(_Cdef_LBFGS_CPP, self).get_result()
        res["H0"] = float(lbfgs_ptr.get_H0())
        return res

class LBFGS_CPP(_Cdef_LBFGS_CPP):
    """This class is the python interface for the c++ LBFGS implementation
    """
