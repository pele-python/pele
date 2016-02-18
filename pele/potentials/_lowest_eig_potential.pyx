"""
# distutils: language = C++
"""
import numpy as np

cdef class cppOrthogonalizeTranslational(object):
    """
    Python interface to C++ implementation of OrthogonalizeTranslational.

    Parameters
    ----------
    natoms : number of atoms

    ndim : integer
        Euclidean dimension of simulation box

    tol : double
        gramm schmidt tolerance
    """
    cdef shared_ptr[cOrthogonalizeTranslational] thisptr
    def __cinit__(self, size_t natoms, size_t ndim, double tol=1e-6):
        self.thisptr = shared_ptr[cOrthogonalizeTranslational](<cOrthogonalizeTranslational*> new cOrthogonalizeTranslational(natoms, ndim, tol))

    def orthogonalize(self, vector):
        """
        :param vector: np.array
            vector to orthogonalize
        :return: void
        """
        self.thisptr.get().orthogonalize(array_wrap_np(np.zeros(vector.size)),
                                         array_wrap_np(vector))

# note that unless I specify the full path pele.potentials._pele.BasePotential
# cython complains about it. Strangely if I move this .pyx and .pxd files in pele/potentials
# the compiler does not complain about incomplete data type
cdef class cppLowestEigPotential(_pele.BasePotential):
    """
    Cython interface to C++ implementation of LowestEigPotential.

    Parameters
    ----------
    potential : _pele.BasePotential
        potential of underlying energy landscape

    coords : np.array
        stationary point for which the lowest eigenvalue is computed

    ndim : integer
        Euclidean dimension of simulation box

    d : double
        finite difference step
    """
    cdef cLowestEigPotential* newptr
    cdef _pele.BasePotential potptr
    def __cinit__(self, landscape_pot, coords, ndim, d=1e-6):
        self.potptr = <_pele.BasePotential> landscape_pot
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cLowestEigPotential(self.potptr.thisptr, array_wrap_np(coords), ndim, d))
        self.newptr = <cLowestEigPotential*> self.thisptr.get()

    def reset_coords(self, coords):
        self.newptr.reset_coords(array_wrap_np(coords))