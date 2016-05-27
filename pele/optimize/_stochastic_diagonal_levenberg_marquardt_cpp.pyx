# distutils: language = c++
from cpython cimport bool as cbool
import numpy as np
cimport numpy as np

from pele.potentials import _pele
from pele.potentials._pythonpotential import as_cpp_potential

cdef class _Cdef_StochasticDiagonalLevenbergMarquardt(_pele_opt.GradientOptimizer):
    """
    Stochastic diagonal Levenberg Marquardt
    
    For a cost function which is a sum of N terms:
    E(x) = \sum_{j=1}^N E_j(x)
    At each iteration, choose one term (say E_k) uniformly at random and do gradient descent:
    (x_{n+1})_i = (x_{n})_i - \eta_{i, n} (\grad E_k(x_{n}))_i
    In addition, also learn the stepsize for each coordinate \eta_{i, n}
    from the second derivatives of the cost function.
    
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
    epsilon : double (optional)
        Global learning rate
    mu : double (optional)
        Parameter to prevent stepsize from blowing up if second
        derivative gets small. (Maximum eta is epsilon/mu.)
    gamma : double (optional)
        Memory parameter for learning of stepsize ('small').
        Larger gamma means less memory.
    tol : double (optional)
        Tolerance on gradient RSM
    nsteps : integer (optional)
        Maximum number of iterations
    seed : integer (optional)
        Random number generator seed
    verbose: bool (optional)
        Print information
    """
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
