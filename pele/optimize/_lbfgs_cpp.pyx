"""
# distutils: language = C++
"""
import sys

import numpy as np

from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
from pele.optimize import Result
from pele.potentials import _pythonpotential

cimport numpy as np
cimport pele.optimize._pele_opt as _pele_opt
from pele.optimize._pele_opt cimport shared_ptr
cimport cython
from cpython cimport bool as cbool

# import the externally defined ljbfgs implementation
cdef extern from "pele/lbfgs.h" namespace "pele":
    cdef cppclass cppLBFGS "pele::LBFGS":
        cppLBFGS(shared_ptr[_pele.cBasePotential], _pele.Array[double], double, int) except +

        void set_H0(double) except +
        void set_tol(double) except +
        void set_maxstep(double) except +
        void set_max_f_rise(double) except +
        void set_use_relative_f(int) except +
        void set_max_iter(int) except +
        void set_iprint(int) except +
        void set_verbosity(int) except +

        double get_H0() except +



cdef class _Cdef_LBFGS_CPP(_pele_opt.GradientOptimizer):
    """This class is the python interface for the c++ LBFGS implementation
    """
    cdef _pele.BasePotential pot
    
    def __cinit__(self, x0, potential, double tol=1e-5, int M=4, double maxstep=0.1, 
                  double maxErise=1e-4, double H0=0.1, int iprint=-1,
                  energy=None, gradient=None,
                  int nsteps=10000, int verbosity=0, events=None, logger=None,
                  rel_energy=False):
        if not issubclass(potential.__class__, _pele.BasePotential):
            if verbosity > 0:
                print "LBFGS_CPP: potential is not subclass of BasePotential; wrapping it.", potential
            potential = _pythonpotential.CppPotentialWrapper(potential)

        self.pot = potential
        if logger is not None:
            print "warning c++ LBFGS is ignoring logger"
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                new cppLBFGS(self.pot.thisptr, 
                             _pele.Array[double](<double*> x0c.data, x0c.size),
                             tol, M) )
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
    
    def get_result(self):
        cdef cppLBFGS* lbfgs_ptr = <cppLBFGS*> self.thisptr.get()
        res = super(_Cdef_LBFGS_CPP, self).get_result()
        res["H0"] = float(lbfgs_ptr.get_H0())
        return res

class LBFGS_CPP(_Cdef_LBFGS_CPP):
    """This class is the python interface for the c++ LBFGS implementation
    """
