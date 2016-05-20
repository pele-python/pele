"""
# distutils: language = C++
"""
import numpy as np
cimport numpy as np
import copy as c
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np_size_t

cdef extern from "pele/xy_model_online.h" namespace "pele":
    cdef cppclass cppXYModelOnline "pele::XYModelOnline":
        cppXYModelOnline(size_t, _pele.Array[size_t], _pele.Array[size_t]) except+
        
cdef class _Cdef_XYModelOnline(_pele.BasePotential):
    """
    Python interface to XYModelOnline.
    
    The potential energy function is:
    U(x) = -\sum_{v\in\text{vertices}} \sum_{e\in\text{edges}(v)} \cos[x(v) - x(v + e)]
    
    Parameters
    ----------
    nr_spins : integer
        Number of spins (vertices) in the graph
        
    links : array of 1d arrays of length 2
        List of links in the graph
    """
    cdef _pele.Array[size_t] head
    cdef _pele.Array[size_t] tail
    def __cinit__(self, nr_spins, links):
        head_ = c.deepcopy(np.array(links[:, 0]))
        tail_ = c.deepcopy(np.array(links[:, 1]))
        head = array_wrap_np_size_t(head_)
        tail = array_wrap_np_size_t(tail_)
        self.thisptr = shared_ptr[_pele.cBasePotential](
            <_pele.cBasePotential*> new cppXYModelOnline(nr_spins, head,
            tail))

class XYModelOnline(_Cdef_XYModelOnline):
    """
    Interface to C++ implementation
    """
