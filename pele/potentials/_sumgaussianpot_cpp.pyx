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

cdef extern from "pele/sumgaussianpot.h" namespace "pele":
    cdef cppclass cppGaussianPot "pele::GaussianPot":
        cppGaussianPot(_pele.Array[double], _pele.Array[double]) except+

cdef class GaussianPot(_pele.BasePotential):
    """python interface to c++ GaussianPot
    """
    def __cinit__(self, means, cov):
        if (means.shape != cov.shape):
            print("means.shape", means.shape)
            print("cov.shape", cov.shape)
            raise Exception("GaussianPot: illegal input")
        means = means.flatten()
        cov = cov.flatten()
        cdef _pele.Array[double] cmeans = array_wrap_np(means)
        cdef _pele.Array[double] ccov = array_wrap_np(cov)
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cppGaussianPot(cmeans, ccov))

cdef class SumGaussianPot(_pele.BasePotential):
    """Sum of GaussianPot
    """
    def __cinit__(self, means, cov):
        if (means.shape != cov.shape):
            print("means.shape", means.shape)
            print("cov.shape", cov.shape)
            raise Exception("SumGaussianPot: illegal input")
        cdef _pele.cCombinedPotential* combpot = new _pele.cCombinedPotential()
        ngauss = means.shape[0]
        bdim = means.shape[1]
        means = means.flatten()
        cov = cov.flatten()
        for i in xrange(ngauss):
            pot = GaussianPot(means[i * bdim : (i + 1) * bdim], cov[i * bdim : (i + 1) * bdim])
            combpot.add_potential(pot.thisptr)
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> combpot)
