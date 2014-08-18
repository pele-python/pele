"""
# distutils: language = C++
"""
import numpy as np
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
cimport numpy as np
from cpython cimport bool
from ctypes import c_size_t as size_t

# cython has no support for integer template argument.  This is a hack to get around it
# https://groups.google.com/forum/#!topic/cython-users/xAZxdCFw6Xs
# Basically you fool cython into thinking INT2 is the type integer,
# but in the generated c++ code you use 2 instead.
# The cython code MyClass[INT2] will create c++ code MyClass<2>.
cdef extern from *:
    ctypedef int INT2 "2"    # a fake type
    ctypedef int INT3 "3"    # a fake type

# use external c++ class
cdef extern from "pele/hs_wca.h" namespace "pele":
    cdef cppclass  cHS_WCA "pele::HS_WCA"[ndim]:
        cHS_WCA(double eps, double sca, _pele.Array[double] radii) except +
    cdef cppclass  cHS_WCAPeriodic "pele::HS_WCAPeriodic"[ndim]:
        cHS_WCAPeriodic(double eps, double sca, _pele.Array[double] radii, _pele.Array[double] boxvec) except +
    cdef cppclass  cHS_WCANeighborList "pele::HS_WCANeighborList":
        cHS_WCANeighborList(_pele.Array[long] & ilist, double eps, double sca, _pele.Array[double] radii) except +    
    cdef cppclass  cHS_WCAFrozen "pele::HS_WCAFrozen"[ndim]:
        cHS_WCAFrozen(double eps, double sca, _pele.Array[double] radii, _pele.Array[double]& reference_coords, _pele.Array[size_t]& frozen_dof) except +
    cdef cppclass  cHS_WCAPeriodicFrozen "pele::HS_WCAPeriodicFrozen"[ndim]:
        cHS_WCAPeriodicFrozen(double eps, double sca, _pele.Array[double] radii, _pele.Array[double] boxvec, _pele.Array[double]& reference_coords, _pele.Array[size_t]& frozen_dof) except +

cdef class HS_WCA(_pele.BasePotential):
    """define the python interface to the c++ HS_WCA implementation
    """
    cpdef bool periodic 
    #thickness of the wca shell is sca * R where R is the hard core radius of the sphere
    def __cinit__(self, eps, sca, radii, ndim=3, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * ndim
        cdef np.ndarray[double, ndim=1] bv
        cdef np.ndarray[double, ndim=1] radiic = np.array(radii, dtype=float)    
        
        if boxvec is None:
            self.periodic = False
            if ndim == 2:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cHS_WCA[INT2](eps, sca, _pele.Array[double](<double*> radiic.data, radiic.size)) )
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cHS_WCA[INT3](eps, sca, _pele.Array[double](<double*> radiic.data, radiic.size)) )
        else:
            self.periodic = True
            ndim = len(boxvec)
            bv = np.array(boxvec, dtype=float)
            if ndim == 2:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                         cHS_WCAPeriodic[INT2](eps, sca, _pele.Array[double](<double*> radiic.data, radiic.size), _pele.Array[double](<double*> bv.data, bv.size)) )
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new 
                         cHS_WCAPeriodic[INT3](eps, sca, _pele.Array[double](<double*> radiic.data, radiic.size), _pele.Array[double](<double*> bv.data, bv.size)) )

cdef class HS_WCAFrozen(_pele.BasePotential):
    """define the python interface to the c++ HS_WCAFrozen implementation
    """
    cpdef bool periodic
    def __cinit__(self, np.ndarray[double, ndim=1] reference_coords, frozen_atoms, eps, sca, np.ndarray[double, ndim=1] radii, ndim=3, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
        cdef np.ndarray[long, ndim=1] frozen_dof
        frozen_dof = np.array([range(ndim*i,ndim*i+ndim) for i in frozen_atoms], dtype=int).reshape(-1) 
                    
        if boxl is not None:
            boxvec = np.array([boxl] * ndim)
        cdef np.ndarray[double, ndim=1] bv
          
        if boxvec is None:
            self.periodic = False
            if ndim==2:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
                     cHS_WCAFrozen[INT2](eps, sca, _pele.Array[double](<double*> radii.data, radii.size), 
                                     _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                                     _pele.Array[size_t](<size_t *> frozen_dof.data, frozen_dof.size) ) )
            elif ndim==3:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
                     cHS_WCAFrozen[INT3](eps, sca, _pele.Array[double](<double*> radii.data, radii.size), 
                                   _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                                   _pele.Array[size_t](<size_t *> frozen_dof.data, frozen_dof.size) ) )
            else:
                raise Exception("HS_WCAFrozen: illegal ndim")
        else:
            self.periodic = True
            bv = np.array(boxvec, dtype=float)
            assert bv.size == ndim
            if ndim==2:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
                     cHS_WCAPeriodicFrozen[INT2](eps, sca, _pele.Array[double](<double*> radii.data, radii.size), 
                                             _pele.Array[double](<double*> bv.data, bv.size), 
                                             _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                                             _pele.Array[size_t](<size_t *> frozen_dof.data, frozen_dof.size) ) )
            elif ndim==3:
                self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
                     cHS_WCAPeriodicFrozen[INT3](eps, sca, _pele.Array[double](<double*> radii.data, radii.size),
                                           _pele.Array[double](<double*> bv.data, bv.size), 
                                           _pele.Array[double](<double *> reference_coords.data, reference_coords.size),
                                           _pele.Array[size_t](<size_t *> frozen_dof.data, frozen_dof.size) ) )
            else:
                raise Exception("HS_WCAFrozen: illegal ndim")

cdef class HS_WCANeighborList(_pele.BasePotential):
    """define the python interface to the c++ HS_WCA implementation
    """
    def __cinit__(self, np.ndarray[long, ndim=1] ilist, eps, sca, np.ndarray[double, ndim=1] radii):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
             cHS_WCANeighborList( _pele.Array[long](<long*> ilist.data, <int> ilist.size), 
                                  eps, sca, 
                                  _pele.Array[double](<double*> radii.data, radii.size)) )

