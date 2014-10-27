"""
# distutils: language = C++
"""
import numpy as np
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from pele.potentials._pele cimport array_wrap_np
from pele.potentials._pele cimport array_wrap_np_size_t
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
    cdef cppclass  cHS_WCACellLists "pele::HS_WCACellLists"[ndim]:
        cHS_WCACellLists(double eps, double sca, _pele.Array[double] radii,
                         _pele.Array[double] boxvec, _pele.Array[double] coords,
                         double rcut, double ncellsx_scale) except +
    cdef cppclass  cHS_WCAPeriodic "pele::HS_WCAPeriodic"[ndim]:
        cHS_WCAPeriodic(double eps, double sca, _pele.Array[double] radii,
                        _pele.Array[double] boxvec) except +
    cdef cppclass  cHS_WCAPeriodicCellLists "pele::HS_WCAPeriodicCellLists"[ndim]:
        cHS_WCAPeriodicCellLists(double eps, double sca,
                                 _pele.Array[double] radii, _pele.Array[double] boxvec, 
                                 _pele.Array[double] coords, double rcut,
                                 double ncellx_scale) except +
    cdef cppclass  cHS_WCANeighborList "pele::HS_WCANeighborList":
        cHS_WCANeighborList(_pele.Array[long] & ilist, double eps, double sca,
                            _pele.Array[double] radii) except +    
    cdef cppclass  cHS_WCAFrozen "pele::HS_WCAFrozen"[ndim]:
        cHS_WCAFrozen(double eps, double sca, _pele.Array[double] radii,
                      _pele.Array[double]& reference_coords,
                      _pele.Array[size_t]& frozen_dof) except +
    cdef cppclass  cHS_WCAPeriodicFrozen "pele::HS_WCAPeriodicFrozen"[ndim]:
        cHS_WCAPeriodicFrozen(double eps, double sca, _pele.Array[double] radii,
                              _pele.Array[double] boxvec,
                              _pele.Array[double]& reference_coords,
                              _pele.Array[size_t]& frozen_dof) except +
    cdef cppclass  cHS_WCACellListsFrozen "pele::HS_WCACellListsFrozen"[ndim]:
        cHS_WCACellListsFrozen(double eps, double sca, _pele.Array[double] radii,
                               _pele.Array[double] boxvec,
                               _pele.Array[double]& reference_coords,
                               _pele.Array[size_t]& frozen_dof, double rcut,
                               double ncellx_scale) except +
    cdef cppclass  cHS_WCAPeriodicCellListsFrozen "pele::HS_WCAPeriodicCellListsFrozen"[ndim]:
        cHS_WCAPeriodicCellListsFrozen(double eps, double sca,
                                       _pele.Array[double] radii,
                                       _pele.Array[double] boxvec,
                                       _pele.Array[double] reference_coords,
                                       _pele.Array[size_t]& frozen_dof,
                                       double rcut, double ncellx_scale)

cdef class HS_WCA(_pele.BasePotential):
    """
    Define the python interface to the c++ HS_WCA implementation.
    
    Should work for most cases, more specialised versions below.
    
    Prameters
    ---------
    eps : float
        wca parameter
    sca : float
        the thickness of the wca shell is sca * R where R is the hard core
        radius of the sphere
    radii : array
        list of radii of the particles
    ndim : integer
        Euclidean dimension of simulation box
    use_periodic : bool, optinal
        Flag to switch on periodic boundary conditions
    boxvec : array
        Box vector
    boxl : float
        In case the box is a cube, the cube length can be given as boxl instead
        of providing boxvec
    use_frozen : bool
        Flag to switch on freezing of selected particles
    reference_coords : array
        Initial values of all coordinates, consisting of frozen and non-frozen
        dof
    frozen_atoms : array
        List of labels of frozen (immobile) particles. Note: This is not the
        list of frozen degrees of freedom
    use_cell_lists : bool
        Flag to switch on usage of cell lists
    rcut : float
        Cutoff radius of cell lists
    ncellx_scale : float
        Parameter controlling the cell list grid spacing: values larger than
        unity lead to finer cell meshing
    """
    cpdef bool periodic
    def __cinit__(self, eps=1.0, sca=1.2,
                  np.ndarray[double, ndim=1] radii=None, ndim=3, boxvec=None,
                  boxl=None, use_periodic=False, use_frozen=False,
                  use_cell_lists=False,
                  np.ndarray[double, ndim=1] reference_coords=None,
                  frozen_atoms=None, rcut=None, ncellx_scale=1.0):
        assert not (boxvec is not None and boxl is not None)
        cdef np.ndarray[long, ndim=1] frozen_dof 
        if boxl is not None:
            boxvec = np.array([boxl] * ndim)
        cdef np.ndarray[double, ndim=1] bv
        if use_frozen:
            if frozen_atoms is None:
                print "HS_WCA: warning: initialising frozen particle potential without frozen particles"
            if reference_coords is None:
                raise Exception("missing input: can not initialise frozen particle potential without reference coordinates")
            else:
                frozen_dof = np.array([range(ndim * i, ndim * i + ndim) for i in frozen_atoms], dtype=int).reshape(-1)
        if use_cell_lists and (rcut is None or boxvec is None):
            raise Exception("HS_WCA: illegal input")
        bv = None
        if use_cell_lists or use_periodic:
            if boxvec is None:
                raise Exception("boxvec is not specified")
            bv = np.array(boxvec, dtype=float)
            assert bv.size == ndim
        cdef _pele.Array[double] rd_
        cdef _pele.Array[double] bv_
        cdef _pele.Array[double] rc_
        cdef _pele.Array[size_t] fd_
        if radii is not None:
            rd_ = array_wrap_np(radii)
        if bv is not None:
            bv_ = array_wrap_np(bv)
        if reference_coords is not None:
            rc_ = array_wrap_np(reference_coords)
        if frozen_atoms is not None:
            fd_ = array_wrap_np_size_t(frozen_dof)
        if use_frozen:
            """
            frozen
            """
            if not use_periodic:
                self.periodic = False
                if ndim == 2:
                    if not use_cell_lists:
                        """
                        frozen, 2d, cartesian, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCAFrozen[INT2](eps, sca, rd_, rc_, fd_))
                    else:
                        """
                        frozen, 2d, cartesian, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCACellListsFrozen[INT2](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
                elif ndim == 3:
                    if not use_cell_lists:
                        """
                        frozen, 3d, cartesian, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCAFrozen[INT3](eps, sca, rd_, rc_, fd_))
                    else:
                        """
                        frozen, 3d, cartesian, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCACellListsFrozen[INT3](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
                else:
                    raise Exception("HS_WCAFrozen: illegal ndim")
            else:
                self.periodic = True
                if ndim == 2:
                    if not use_cell_lists:
                        """
                        frozen, 2d, periodic, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCAPeriodicFrozen[INT2](eps, sca, rd_, bv_, rc_, fd_))
                    else:
                        """
                        frozen, 2d, periodic, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                       cHS_WCAPeriodicCellListsFrozen[INT2](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
                elif ndim == 3:
                    if not use_cell_lists:
                        """
                        frozen, 3d, periodic, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCAPeriodicFrozen[INT3](eps, sca, rd_, bv_, rc_, fd_))
                    else:
                        """
                        frozen, 3d, periodic, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                       cHS_WCAPeriodicCellListsFrozen[INT3](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
                else:
                    raise Exception("HS_WCAFrozen: illegal ndim")
        else:
            """
            non-frozen
            """
            assert(not use_frozen)
            if not use_periodic:
                self.periodic = False
                if ndim == 2:
                    if not use_cell_lists:
                        """
                        non-frozen, 2d, cartesian, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCA[INT2](eps, sca, rd_))
                    else:
                        """
                        non-frozen, 2d, cartesian, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                             cHS_WCACellLists[INT2](eps, sca, rd_, bv_, rc_, rcut, ncellx_scale))
                elif ndim == 3:
                    if not use_cell_lists:
                        """
                        non-frozen, 3d, cartesian, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                             cHS_WCA[INT3](eps, sca, rd_))
                    else:
                        """
                        non-frozen, 3d, cartesian, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCACellLists[INT3](eps, sca, rd_, bv_, rc_, rcut, ncellx_scale))
                else:
                    raise Exception("HS_WCA: illegal ndim")
            else:
                self.periodic = True
                if ndim == 2:
                    if not use_cell_lists:
                        """
                        non-frozen, 2d, periodic, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                             cHS_WCAPeriodic[INT2](eps, sca, rd_, bv_))
                    else:
                        """
                        non-frozen, 2d, periodic, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                             cHS_WCAPeriodicCellLists[INT2](eps, sca, rd_, bv_, rc_, rcut, ncellx_scale))
                elif ndim == 3:
                    if not use_cell_lists:
                        """
                        non-frozen, 3d, periodic, no cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                             cHS_WCAPeriodic[INT3](eps, sca, rd_, bv_))
                    else:
                        """
                        non-frozen, 3d, periodic, use cell lists
                        """
                        self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                             cHS_WCAPeriodicCellLists[INT3](eps, sca, rd_, bv_, rc_, rcut, ncellx_scale)) 
                else:
                    raise Exception("HS_WCA: illegal ndim")
            

cdef class HS_WCASimple(_pele.BasePotential):
    """define the python interface to the c++ HS_WCA implementation
    
    Parameters
    ----------
    eps : float
        wca parameter
    sca : float
        the thickness of the wca shell is sca * R where R is the hard core radius of the sphere
    radii : array
        list of radii of the particles
    """

    cpdef bool periodic 
    #thickness of the wca shell is sca * R where R is the hard core radius of the sphere
    def __cinit__(self, eps, sca, radii, ndim=3, boxvec=None, boxl=None):
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * ndim
        cdef np.ndarray[double, ndim=1] bv
        cdef np.ndarray[double, ndim=1] radiic = np.array(radii, dtype=float)
        rd_ = array_wrap_np(radiic)    
        if boxvec is None:
            self.periodic = False
            if ndim == 2:
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*>new
                    cHS_WCA[INT2](eps, sca, rd_))
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*>new
                    cHS_WCA[INT3](eps, sca, rd_))
        else:
            self.periodic = True
            ndim = len(boxvec)
            bv = np.array(boxvec, dtype=float)
            bv_ = array_wrap_np(bv)
            if ndim == 2:
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*>new 
                         cHS_WCAPeriodic[INT2](eps, sca, rd_, bv_))
            else:
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*>new 
                         cHS_WCAPeriodic[INT3](eps, sca, rd_, bv_))

cdef class HS_WCAFrozen(_pele.BasePotential):
    """define the python interface to the c++ HS_WCAFrozen implementation
    """
    cpdef bool periodic
    def __cinit__(self, np.ndarray[double, ndim=1] reference_coords,
                  frozen_atoms, np.ndarray[double, ndim=1] radii, eps=1.0,
                  sca=1.2, ndim=3, use_periodic=False, boxvec=None, boxl=None,
                  use_cell_lists=False, rcut=None, ncellx_scale=1.0):
        assert not (boxvec is not None and boxl is not None)
        cdef np.ndarray[long, ndim=1] frozen_dof
        frozen_dof = np.array([range(ndim*i,ndim*i+ndim) for i in frozen_atoms], dtype=int).reshape(-1) 
                    
        if boxl is not None:
            boxvec = np.array([boxl] * ndim)
        cdef np.ndarray[double, ndim=1] bv
        
        if use_cell_lists and (rcut is None or boxvec is None):
            raise Exception("HS_WCAFrozen: illegal input")
          
        bv = None
        if use_cell_lists or use_periodic:
            if boxvec is None:
                raise Exception("boxvec is not specified")
            bv = np.array(boxvec, dtype=float)
            assert bv.size == ndim
        
        rd_ = array_wrap_np(radii)
        rc_ = array_wrap_np(reference_coords)
        fd_ = array_wrap_np_size_t(frozen_dof)
        cdef _pele.Array[double] bv_
        if bv is not None:
            bv_ = array_wrap_np(bv)
        if not use_periodic:
            self.periodic = False
            if ndim==2:
                if not use_cell_lists:
                    """
                    frozen, 2d, cartesian, no cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                         cHS_WCAFrozen[INT2](eps, sca, rd_, rc_, fd_))
                else:
                    """
                    frozen, 2d, cartesian, use cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                         cHS_WCACellListsFrozen[INT2](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
            elif ndim==3:
                if not use_cell_lists:
                    """
                    frozen, 3d, cartesian, no cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                         cHS_WCAFrozen[INT3](eps, sca, rd_, rc_, fd_))
                else:
                    """
                    frozen, 3d, cartesian, use cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                         cHS_WCACellListsFrozen[INT3](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
            else:
                raise Exception("HS_WCAFrozen: illegal ndim")
        else:
            self.periodic = True
            if ndim==2:
                if not use_cell_lists:
                    """
                    frozen, 2d, periodic, no cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                         cHS_WCAPeriodicFrozen[INT2](eps, sca, rd_, bv_, rc_, fd_))
                else:
                    """
                    frozen, 2d, periodic, use cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                   cHS_WCAPeriodicCellListsFrozen[INT2](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
            elif ndim==3:
                if not use_cell_lists:
                    """
                    frozen, 3d, periodic, no cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new 
                         cHS_WCAPeriodicFrozen[INT3](eps, sca, rd_, bv_, rc_, fd_))
                else:
                    """
                    frozen, 3d, periodic, use cell lists
                    """
                    self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new
                                   cHS_WCAPeriodicCellListsFrozen[INT3](eps, sca, rd_, bv_, rc_, fd_, rcut, ncellx_scale))
            else:
                raise Exception("HS_WCAFrozen: illegal ndim")

cdef class HS_WCAPeriodicCellLists(_pele.BasePotential):
    """define the python interface to the c++ HS_WCAPeriodicCellLists implementation
    """
    cpdef bool frozen
    def __cinit__(self, eps, sca, radii, boxvec, coords, rcut, ndim=3,
                  ncellx_scale=1.0, frozen_atoms=None):
            ndim = len(boxvec)
            cdef np.ndarray[double, ndim=1] radiic = np.array(radii, dtype=float)
            rd_ = array_wrap_np(radiic)
            cdef np.ndarray[double, ndim=1] boxvecc = np.array(boxvec, dtype=float)
            cdef np.ndarray[double, ndim=1] coordsc = np.array(coords, dtype=float)
            cdef np.ndarray[long, ndim=1] frozen_dof
            bv_ = array_wrap_np(boxvecc)
            co_ = array_wrap_np(coordsc)
            if frozen_atoms is None:
                self.frozen = False
                if ndim == 2:
                    self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new
                                   cHS_WCAPeriodicCellLists[INT2](eps, sca, rd_, bv_, co_, rcut, ncellx_scale))
                elif ndim == 3:
                    self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new
                                   cHS_WCAPeriodicCellLists[INT3](eps, sca, rd_, bv_, co_, rcut, ncellx_scale))
                else:
                    raise Exception("HS_WCAPeriodicCellLists: illegal boxdimension")
            else:
                self.frozen = True
                frozen_dof = np.array([range(ndim * i, ndim * i + ndim) for i in frozen_atoms], dtype = int).reshape(-1) 
                fd_ = array_wrap_np_size_t(frozen_dof)
                if ndim == 2:
                    self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new
                                   cHS_WCAPeriodicCellListsFrozen[INT2](eps, sca, rd_, bv_, co_, fd_, rcut, ncellx_scale)) 
                elif ndim == 3:
                    self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new
                                   cHS_WCAPeriodicCellListsFrozen[INT3](eps, sca, rd_, bv_, co_, fd_, rcut, ncellx_scale))
                else:
                    raise Exception("HS_WCAPeriodicCellLists: illegal boxdimension")

cdef class HS_WCANeighborList(_pele.BasePotential):
    """define the python interface to the c++ HS_WCA implementation
    """
    def __cinit__(self, np.ndarray[long, ndim=1] ilist, eps, sca, np.ndarray[double, ndim=1] radii):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*> new 
             cHS_WCANeighborList(_pele.Array[long](<long*> ilist.data, <int> ilist.size),
                                 eps, sca, array_wrap_np(radii)))
