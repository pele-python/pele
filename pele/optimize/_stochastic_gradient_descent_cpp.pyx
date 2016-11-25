#distutils: language = c++
from cpython cimport bool as cbool
import numpy as np
cimport numpy as np

from pele.potentials import _pele
from pele.potentials._pythonpotential import as_cpp_potential

cdef class _Cdef_StochasticGradientDescent(_pele_opt.GradientOptimizer):
    """
    Stochastic gradient descent
    
    For a cost function which is a sum of N terms:
    E(x) = \sum_{j=1}^N E_j(x)
    At each iteration, choose one term (say E_k) uniformly at random and do gradient descent:
    x_{n+1} = x_{n} - \eta \grad E_k(x_{n})
    
    References
    ----------
    * Christopher M. Bishop: Pattern Recognition and Machine Learning
    * Yann LeCun et al.: http://yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf

    
    Parameters
    ----------
    x0 : array of doubles
        Initial coordinates
    potential : pele potential
        Cost function
    eta : double (optional)
        Stepsize
    nsteps : integer (optional)
        Maximum number of iterations
    tol : double (optional)
        Tolerance on gradient RSM
    seed : integer (optional)
        Random number generator seed
    """
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
