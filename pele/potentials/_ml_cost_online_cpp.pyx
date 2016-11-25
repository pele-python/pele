"""
# distutils: language = C++
"""
import numpy as np
cimport numpy as np
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np

cdef extern from "pele/ml_cost_online.h" namespace "pele":
    cdef cppclass cppMLCostOnlineGauss "pele::MLCostOnlineGauss":
        cppMLCostOnlineGauss(_pele.Array[double]) except+
        
cdef class _Cdef_MLCostOnline(_pele.BasePotential):
    """
    Python interface to MLCostOnline.
    
    Parameters
    ----------
    data : array of doubles
        Recorded data points.
    
    model : string (optional)
        Name of probability distribution function to fit to the data
    """
    cdef _pele.Array[double] datac
    def __cinit__(self, data, model="gauss"):
        datac = array_wrap_np(data)
        if model == "gauss":
            self.thisptr = shared_ptr[_pele.cBasePotential](
                <_pele.cBasePotential*> new cppMLCostOnlineGauss(datac))
        else:
            raise Exception("MLCostOnline: illegal input: model name")

class MLCostOnline(_Cdef_MLCostOnline):
    """
    Interface to C++ implementation
    """
