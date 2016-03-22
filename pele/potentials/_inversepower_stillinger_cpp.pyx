"""
# distutils: language = C++
"""
import numpy as np

cimport numpy as np
from cpython cimport bool

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np

# https://groups.google.com/forum/#!topic/cython-users/xAZxdCFw6Xs
cdef extern from *:
    ctypedef int INT2 "2"    # a fake type
    ctypedef int INT3 "3"    # a fake type

cdef extern from "pele/inversepower_stillinger.h" namespace "pele":
    cdef cppclass cInversePowerStillinger "pele::InversePowerStillinger"[ndim]:
        cInversePowerStillinger(size_t pow, _pele.Array[double] radii) except +
    cdef cppclass cInversePowerStillingerPeriodic "pele::InversePowerStillingerPeriodic"[ndim]:
        cInversePowerStillingerPeriodic(size_t pow, _pele.Array[double] radii, _pele.Array[double] boxvec) except +
        
cdef class InversePowerStillinger(_pele.BasePotential):
    """
    Python interface to C++ implementation of InversePowerStillinger.
    
    Parameters
    ----------
    pow : integer
        Exponent value
        
    radii : np.array
        List of particles radii
    
    ndim : integer
        Euclidean dimension of simulation box
    
    boxvec : array
        Box vector
    
    boxl : float
        In case the box is a cube, the cube length can be given as boxl
        instead of providing boxvec
    """
    cpdef bool periodic
    cdef _pele.Array[double] bv_, radii_
    def __cinit__(self, pow, radii, ndim=3, boxvec=None, boxl=None):
        assert(ndim == 2 or ndim == 3)
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * ndim
        radii_ = array_wrap_np(radii)
        if boxvec is not None:
            if len(boxvec) != ndim:
                raise Exception("InversePowerStillinger: illegal input, illegal boxvec")
            bv_ = array_wrap_np(boxvec)
            if ndim == 2:
                # no cell lists, periodic, 2d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillingerPeriodic[INT2](pow, radii_, bv_))
            else:
                # no cell lists, periodic, 3d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillingerPeriodic[INT3](pow, radii_, bv_))
        else:
            if ndim == 2:
                # no cell lists, non-periodic, 2d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillinger[INT2](pow, radii_))
            else:
                # no cell lists, non-periodic, 3d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillinger[INT3](pow, radii_))
