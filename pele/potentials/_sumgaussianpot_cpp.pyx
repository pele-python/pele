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
    def __cinit__(self, np.ndarray[double, ndim=2] means, np.ndarray[double, ndim=2] cov):
        if (means.size() != cov.size()):
            raise Exception("SumGaussianPot: illegal input")
        bdim = means.size()[0]
        cdef _pele.Array[double] m_ = array_wrap_np(np.flatten(means))
        cdef _pele.Array[double] c_ = array_wrap_np(np.flatten(cov))
        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cppSumGaussianPot(bdim, m_, c_))
