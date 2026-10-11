"""
# distutils: language = C++
# cython: language_level=3str
"""
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
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
    ctypedef int INT12 "12"  # a fake type

cdef extern from "pele/inversepower_stillinger_cut_quad.hpp" namespace "pele":
    cdef cppclass cInversePowerStillingerCutQuad "pele::InversePowerStillingerCutQuad"[ndim]:
        cInversePowerStillingerCutQuad(int power, double v0, double cutoff_factor, _pele.Array[double] radii, double non_additivity) except +
    cdef cppclass cInversePowerStillingerCutQuadPeriodic "pele::InversePowerStillingerCutQuadPeriodic"[ndim]:
        cInversePowerStillingerCutQuadPeriodic(int power, double v0, double cutoff_factor, _pele.Array[double] radii, _pele.Array[double] boxvec, double non_additivity) except +
    cdef cppclass cInversePowerStillingerCutQuadPeriodicCellLists "pele::InversePowerStillingerCutQuadPeriodicCellLists"[ndim]:
        cInversePowerStillingerCutQuadPeriodicCellLists(int power, double v0, double cutoff_factor, _pele.Array[double] radii, _pele.Array[double] boxvec, double ncellx_scale, double non_additivity) except +
    cdef cppclass cInversePowerStillingerCutQuadInt "pele::InversePowerStillingerCutQuadInt"[ndim, power]:
        cInversePowerStillingerCutQuadInt(double v0, double cutoff_factor, _pele.Array[double] radii, double non_additivity) except +
    cdef cppclass cInversePowerStillingerCutQuadIntPeriodic "pele::InversePowerStillingerCutQuadIntPeriodic"[ndim, power]:
        cInversePowerStillingerCutQuadIntPeriodic(double v0, double cutoff_factor, _pele.Array[double] radii, _pele.Array[double] boxvec, double non_additivity) except +
    cdef cppclass cInversePowerStillingerCutQuadIntPeriodicCellLists "pele::InversePowerStillingerCutQuadIntPeriodicCellLists"[ndim, power]:
        cInversePowerStillingerCutQuadIntPeriodicCellLists(double v0, double cutoff_factor, _pele.Array[double] radii, _pele.Array[double] boxvec, double ncellx_scale, double non_additivity) except +    


cdef class InversePowerStillingerCutQuad(_pele.BasePotential):
    """
    Python interface to C++ implementation of InversePowerStillingerCutQuad,
    a smooth repulsive potential with a cutoff distance. This corresponds to
    the soft sphere potential in
    https://journals.aps.org/prx/abstract/10.1103/PhysRevX.12.021001.
    For the canonical case of power = 12, this potential delegates to a
    templated C++ implementation for speed.

    Parameters
    ----------
    power : integer
        Exponent value

    v0 : float
        Energy scale.

    cutoff_factor : float
        Cutoff distance as a multiple of the pair interaction length.

    radii : np.array
        List of particles radii

    ndim : integer
        Euclidean dimension of simulation box

    boxvec : array
        Box vector

    boxl : float
        In case the box is a cube, the cube length can be given as boxl
        instead of providing boxvec

    ncellx_scale : float, optional
        Cell-list grid scale.

    use_cell_lists : bool, optional
        Enable cell lists for periodic systems.

    non_additivity : float, optional
        Non-additivity parameter passed to the pair interaction.
    """
    cdef public bool periodic
    cdef _pele.Array[double] bv_, radii_
    def __cinit__(self, power, v0, cutoff_factor, radii, ndim=3, boxvec=None, boxl=None,
                  ncellx_scale=1., use_cell_lists=False, non_additivity=0.):
        assert(ndim == 2 or ndim == 3)
        power = int(power)
        assert not (boxvec is not None and boxl is not None)
        if boxl is not None:
            boxvec = [boxl] * ndim
        boxvec = np.array(boxvec, dtype=np.float64)
        radii_ = array_wrap_np(radii)
        if boxvec is not None:
            if len(boxvec) != ndim:
                raise Exception("InversePowerStillinger: illegal input, illegal boxvec")
            bv_ = array_wrap_np(boxvec)
            if use_cell_lists:
                if ndim == 2 and power == 12:
                    # no cell lists, periodic, 2d
                    self.thisptr = shared_ptr[_pele.cBasePotential](
                        <_pele.cBasePotential*> new cInversePowerStillingerCutQuadIntPeriodicCellLists[INT2, INT12](
                            v0, cutoff_factor, radii_, bv_, ncellx_scale, non_additivity
                            )
                        )
                elif ndim == 3 and power == 12:
                    # no cell lists, periodic, 3d
                    self.thisptr = shared_ptr[_pele.cBasePotential](
                        <_pele.cBasePotential*> new cInversePowerStillingerCutQuadIntPeriodicCellLists[INT3, INT12](
                            v0, cutoff_factor, radii_, bv_, ncellx_scale, non_additivity
                            )
                        )
                elif ndim == 2:
                    self.thisptr = shared_ptr[_pele.cBasePotential](
                        <_pele.cBasePotential*> new cInversePowerStillingerCutQuadPeriodicCellLists[INT2](
                            power, v0, cutoff_factor, radii_, bv_, ncellx_scale, non_additivity
                            ))
                elif ndim == 3:
                    # no cell lists, periodic, 3d
                    self.thisptr = shared_ptr[_pele.cBasePotential](
                        <_pele.cBasePotential*> new cInversePowerStillingerCutQuadPeriodicCellLists[INT3](
                            power, v0, cutoff_factor, radii_, bv_, ncellx_scale, non_additivity
                            )
                        )
                else:
                    raise Exception("InversePowerStillinger: illegal input, illegal ndim")
            else:
                if ndim == 2 and power == 12:
                    # no cell lists, periodic, 2d
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                        cInversePowerStillingerCutQuadIntPeriodic[INT2, INT12](v0, cutoff_factor, radii_, bv_, non_additivity))
                elif ndim == 3 and power == 12:
                    # no cell lists, periodic, 3d
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                        cInversePowerStillingerCutQuadIntPeriodic[INT3, INT12](v0, cutoff_factor, radii_, bv_, non_additivity))
                elif ndim == 2:
                    # no cell lists, periodic, 2d
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                        cInversePowerStillingerCutQuadPeriodic[INT2](power, v0, cutoff_factor, radii_, bv_, non_additivity))
                elif ndim == 3:
                    # no cell lists, periodic, 3d
                    self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                        cInversePowerStillingerCutQuadPeriodic[INT3](power, v0, cutoff_factor, radii_, bv_, non_additivity))
                else:
                    raise Exception("InversePowerStillinger: illegal input, illegal ndim")
        else:
            assert use_cell_lists is False, "Cell Lists Implementation only works for periodic systems"
            if ndim == 2 and power == 12:
                # no cell lists, non-periodic, 2d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillingerCutQuadInt[INT2, INT12](v0, cutoff_factor, radii_, non_additivity))
            elif ndim ==3 and power == 12:
                # no cell lists, non-periodic, 3d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillingerCutQuadInt[INT3, INT12](v0, cutoff_factor, radii_, non_additivity))
            elif ndim == 2:
                # no cell lists, non-periodic, 2d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillingerCutQuad[INT2](power, v0, cutoff_factor, radii_, non_additivity))
            elif ndim ==3:
                # no cell lists, non-periodic, 3d
                self.thisptr = shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new
                                                                                    cInversePowerStillingerCutQuad[INT3](power, v0, cutoff_factor, radii_, non_additivity))
            else:
                raise Exception("InversePowerStillinger: illegal input, illegal ndim")
