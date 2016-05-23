#distutils: language = c++
from cpython cimport bool as cbool
import numpy as np
cimport numpy as np

from pele.potentials import _pele
from pele.potentials._pythonpotential import as_cpp_potential

cdef class _Cdef_StochasticGradientDescent(_pele_opt.GradientOptimizer):
    cdef _pele.BasePotential pot
    cdef cppStochasticGradientDescent* newptr
    def __cinit__(self, x0, potential, eta=0.5, nsteps=10000, tol=1e-5, seed=42,
        verbose=False):
        potential = as_cpp_potential(potential)
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.pot = pot
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer](<_pele_opt.cGradientOptimizer*>
            new cppStochasticGradientDescent(self.pot.thisptr,
                _pele.Array[double](<double*> x0c.data, x0c.size),
                eta, tol, seed, verbose))
        self.newptr = <cppStochasticGradientDescent*> self.thisptr.get()
        self.thisptr.get().set_max_iter(nsteps)

class StochasticGradientDescent(_Cdef_StochasticGradientDescent):
    """
    Interface to C++ implementation.
    """
