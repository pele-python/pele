"""
# distutils: language = C++
"""
import numpy as np
cimport numpy as np
from ctypes import c_size_t as size_t
from ctypes import c_double as double

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport array_wrap_np
from pele.potentials._pele cimport Array
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport BasePotential

cdef extern from "pele/gaussianpot.h" namespace "pele":
    cdef cppclass cppGaussianPot "pele::GaussianPot":
        cppGaussianPot(_pele.Array[double], _pele.Array[double]) except+

cdef class GaussianPot(_pele.BasePotential):
    """python interface to c++ GaussianPot
    
    Gaussian-shaped potential function. For given arrays mean and cov, the energy of a configuration with coordinates $x$ and $N$ degrees of freedom is given by:
    $\text{energy}(\text{mean}; \text{cov}, x) = -\exp \left[ -1/2 \sum_{i=1}^N (x_i - \text{mean}_i)^2 / \text{cov}_i \right]$.
    
    Parameters
    ----------
    mean : array
        Mean of the gaussian, needs to be of lenght bdim.
    cov : array
        Diagonal of the covariance matrix of the gaussian, needs to be
        of length bdim.
    """
    def __cinit__(self, mean, cov):
        if (mean.shape != cov.shape):
            print("mean.shape", mean.shape)
            print("cov.shape", cov.shape)
            raise Exception("GaussianPot: illegal input")
        cdef _pele.Array[double] cmean = array_wrap_np(mean)
        cdef _pele.Array[double] ccov = array_wrap_np(cov)
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cppGaussianPot(cmean, ccov))

cdef class SumGaussianPot(_pele.BasePotential):
    """Sum of GaussianPot
    
    Sum of gaussian-shaped potentials.
    
    Parameters
    ----------
    means : array, (number of gaussians) * (number of degrees of
        freedom)
        Array of means for gaussian potential sum.
    cov : array, (number of gaussians) * (number of degrees of freedom)
        Array of covariance matrix diagonals for gaussian potential sum.
    """
    # The following line is based on the discussion here:
    # https://groups.google.com/forum/#!topic/cython-users/OLrQ0QWQko0
    cdef public bdim
    def __cinit__(self, means, cov):
        if (means.shape != cov.shape):
            print("means.shape", means.shape)
            print("cov.shape", cov.shape)
            raise Exception("SumGaussianPot: illegal input")
        cdef _pele.cCombinedPotential* combpot = new _pele.cCombinedPotential()
        ngauss = means.shape[0]
        self.bdim = means.shape[1]
        for i in xrange(ngauss):
            pot = GaussianPot(means[i,:], cov[i,:])
            combpot.add_potential(pot.thisptr)
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> combpot)
    
    @property
    def dim(self):
        return self.bdim
        
    def getEnergy2ndGradient(self, coords):
        e, g, h = self.getEnergyGradientHessian(coords)
        return e, g, np.diag(h)
