# distutils: language = c++
from cpython cimport bool as cbool
import numpy as np
cimport numpy as np

from pele.potentials import _pele
from pele.potentials._pythonpotential import as_cpp_potential

cdef class _Cdef_SteepestDescent(_pele_opt.GradientOptimizer):
    cdef _pele.BasePotential pot
    cdef cppSteepestDescent* newptr
    def __cinit__(self, x0, potential, eta=1, eta_min=None,
        double tol=1e-4, size_t nsteps=10000, cbool verbose=False):
        potential = as_cpp_potential(potential)
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.pot = pot
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer](<_pele_opt.cGradientOptimizer*>
            new cppSteepestDescent(self.pot.thisptr, _pele.Array[double](<double*> x0c.data, x0c.size),
            eta, tol, verbose))
        self.newptr = <cppSteepestDescent*> self.thisptr.get()
        self.thisptr.get().set_max_iter(nsteps)
        if eta_min is not None:
            self.newptr.set_eta_min(eta_min)

class SteepestDescentCPP(_Cdef_SteepestDescent):
    """
    Interface to C++ implementation.
    """
