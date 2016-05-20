# distutils: language = c++
from cpython cimport bool as cbool
import numpy as np
cimport numpy as np

from pele.potentials import _pele
from pele.potentials._pythonpotential import as_cpp_potential

"""
StochasticDiagonalLevenbergMarquardt(std::shared_ptr<BasePotential> potential,
        const Array<double>& x0, const double epsilon=1,
        const double mu=1, const double gamma=0.1,
        const double tol=1e-5, const size_t seed=42,
        const bool verbose=false)
"""

cdef class _Cdef_StochasticDiagonalLevenbergMarquardt(_pele_opt.GradientOptimizer):
    cdef _pele.BasePotential pot
    cdef cppStochasticDiagonalLevenbergMarquardt* newptr
    def __cinit__(self, x0, potential, epsilon=1, mu=2, gamma=0.1,
        tol=1e-5, nsteps=10000, seed=42, verbose=False):
        potential = as_cpp_potential(potential)
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.pot = pot
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer](<_pele_opt.cGradientOptimizer*>
            new cppStochasticDiagonalLevenbergMarquardt(self.pot.thisptr,
            _pele.Array[double](<double*> x0c.data, x0c.size),
            epsilon, mu, gamma, tol, seed, verbose))
        self.newptr = <cppStochasticDiagonalLevenbergMarquardt*> self.thisptr.get()
        self.thisptr.get().set_max_iter(nsteps)

class StochasticDiagonalLevenbergMarquardt(_Cdef_StochasticDiagonalLevenbergMarquardt):
    """
    Interface to C++ implementation.
    """
