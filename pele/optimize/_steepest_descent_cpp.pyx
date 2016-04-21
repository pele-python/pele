# distutils: language = c++
# distutils: sources = modified_fire.cpp
from cpython cimport bool as cbool
import numpy as np
cimport numpy as np

from pele.potentials import _pele
from pele.potentials._pythonpotential import as_cpp_potential

cdef class _Cdef_SteepestDescent(_pele_opt.GradientOptimizer):
    cdef _pele.BasePotential pot
    def __cinit__(self, x0, potential, double eta=1e-1, double tol=1e-4, size_t nsteps=10000, cbool verbose=False):
        potential = as_cpp_potential(potential)
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer](<_pele_opt.cGradientOptimizer*>
            new cppSteepestDescent(pot.thisptr, _pele.Array[double](<double*> x0c.data, x0c.size),
            eta, tol, verbose))
        self.thisptr.get().set_max_iter(nsteps)
        self.pot = pot

class SteepestDescentCPP(_Cdef_SteepestDescent):
    """
    Interface to C++ implementation.
    """
