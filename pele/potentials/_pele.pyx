"""
# distutils: language = C++

basic potential interface stuff
"""
import numpy as np
cimport numpy as np
from ctypes import c_size_t as size_t


cdef class BasePotential(object):
    """this class defines the python interface for c++ potentials

    Notes
    -----
    for direct access to the underlying c++ potential use self.thisptr
    """

    def getEnergyGradient(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x.size)
        e = self.thisptr.get().get_energy_gradient(array_wrap_np(x),
                                                   array_wrap_np(grad))
        return e, grad

    def getEnergy(self, np.ndarray[double, ndim=1] x not None):
        # redirect the call to the c++ class
        return self.thisptr.get().get_energy(array_wrap_np(x))

    def getGradient(self, np.ndarray[double, ndim=1] x not None):
        e, grad = self.getEnergyGradient(x)
        return grad

    def getEnergyGradientHessian(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] grad = np.zeros(x.size)
        cdef np.ndarray[double, ndim=1] hess = np.zeros(x.size**2)
        e = self.thisptr.get().get_energy_gradient_hessian(array_wrap_np(x),
                                                           array_wrap_np(grad),
                                                           array_wrap_np(hess))
        return e, grad, hess.reshape([x.size, x.size])

    def getHessian(self, np.ndarray[double, ndim=1] x not None):
        cdef np.ndarray[double, ndim=1] hess = np.zeros(x.size**2)
        self.thisptr.get().get_hessian(array_wrap_np(x), array_wrap_np(hess))
        return np.reshape(hess, [x.size, x.size])

    def NumericalDerivative(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] grad = np.zeros([x.size])
        self.thisptr.get().numerical_gradient(array_wrap_np(x),
                                              array_wrap_np(grad),
                                              eps)
        return grad

    def NumericalHessian(self, np.ndarray[double, ndim=1] x not None, double eps=1e-6):
        # redirect the call to the c++ class
        cdef np.ndarray[double, ndim=1] hess = np.zeros([x.size**2])
        self.thisptr.get().numerical_hessian(array_wrap_np(x),
                                             array_wrap_np(hess),
                                             eps)
        return np.reshape(hess, [x.size, x.size])

cdef class PairwisePotentialInterface(BasePotential):
    """this class defines the python interface for c++ pair potentials

    Notes
    -----
    for direct access to the underlying c++ potential use self.thisptr
    """

    def getInteractionGradient(self, double r, int atomi, int atomj):
        cdef double r2 = r * r
        cdef double grad
        (<cPairwisePotentialInterface*>self.thisptr.get()).get_interaction_energy_gradient(r2, &grad, atomi, atomj)
        return grad

    def getInteractionHessian(self, double r, int atomi, int atomj):
        cdef double r2 = r * r
        cdef double grad
        cdef double hess
        (<cPairwisePotentialInterface*>self.thisptr.get()).get_interaction_energy_gradient_hessian(r2, &grad, &hess, atomi, atomj)
        return hess

    def getNeighbors(self, np.ndarray[double, ndim=1] coords not None,
                      include_atoms=None, cutoff_factor=1.0):
        cdef Array[vector[size_t]] c_neighbor_indss
        cdef Array[vector[vector[double]]] c_neighbor_distss
        cdef Array[short] c_include_atoms

        if include_atoms is None:
            (<cPairwisePotentialInterface*>self.thisptr.get()).get_neighbors(
                array_wrap_np(coords), c_neighbor_indss, c_neighbor_distss,
                cutoff_factor)
        else:
            c_include_atoms = Array[short](len(include_atoms))
            for i in xrange(len(include_atoms)):
                c_include_atoms[i] = include_atoms[i]
            (<cPairwisePotentialInterface*>self.thisptr.get()).get_neighbors_picky(
                array_wrap_np(coords), c_neighbor_indss, c_neighbor_distss,
                c_include_atoms, cutoff_factor)

        neighbor_indss = []
        for i in xrange(c_neighbor_indss.size()):
            neighbor_indss.append([])
            for c_neighbor_ind in c_neighbor_indss[i]:
                neighbor_indss[-1].append(c_neighbor_ind)

        neighbor_distss = []
        for i in xrange(c_neighbor_distss.size()):
            neighbor_distss.append([])
            for c_nneighbor_dist in c_neighbor_distss[i]:
                neighbor_distss[-1].append([])
                for dist_comp in c_nneighbor_dist:
                    neighbor_distss[-1][-1].append(dist_comp)

        return neighbor_indss, neighbor_distss

    def getOverlaps(self, np.ndarray[double, ndim=1] coords not None):
        cdef vector[size_t] c_overlap_inds

        c_overlap_inds = (<cPairwisePotentialInterface*>self.thisptr.get()).get_overlaps(
            array_wrap_np(coords))

        overlap_inds = []
        for i in xrange(0, len(c_overlap_inds), 2):
            overlap_inds.append((c_overlap_inds[i], c_overlap_inds[i + 1]))

        return overlap_inds
