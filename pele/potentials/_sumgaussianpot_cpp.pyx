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
    cdef cppclass cppSumGaussianPot "pele::SumGaussianPot":
        cppSumGaussianPot(size_t, _pele.Array[double], _pele.Array[double]) except+

cdef class SumGaussianPot(_pele.BasePotential):
    """python interface to c++ SumGaussianPot
    """
    def __cinit__(self, means, cov):
        if (means.shape != cov.shape):
            print("means.shape", means.shape)
            print("cov.shape", cov.shape)
            raise Exception("SumGaussianPot: illegal input")
        bdim = means.shape[0]
        cdef _pele.Array[double] m_ = array_wrap_np(means.flatten())
        cdef _pele.Array[double] c_ = array_wrap_np(cov.flatten())
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cppSumGaussianPot(bdim, m_, c_))
