"""tools for finding the smallest eigenvalue and associated eigenvector
using Rayleigh-Ritz minimization
"""
from __future__ import print_function

import numpy as np
import logging

from pele.transition_states import orthogopt
from pele.potentials.potential import BasePotential
from pele.optimize import MYLBFGS
import pele.utils.rotations as rotations

__all__ = ["findLowestEigenVector", "analyticalLowestEigenvalue", "FindLowestEigenVector"]


class LowestEigPot(BasePotential):
    """Potential wrapper for use in an optimizer for finding the eigenvector

    Notes
    -----
    here the energy corresponds to the eigenvalue, and the coordinates to be
    optimized is the eigenvector

    Parameters
    -----------
    coords : array
        The point in space where we want to compute the lowest eigenvector
    pot : potential object
        The potential of the system.  i.e. pot.getEnergyGradient(coords)
        gives the energy and gradient
    
    orthogZeroEigs : callable
        The function which makes a vector orthogonal to the known
        eigenvectors with zero eigenvalues.  The default assumes global
        translational and rotational symmetry
    dx : float
        the local curvature is approximated using points separated by dx
    first_order : bool
        use the first order forward finite differences approximation for
        the curvature rather than the second order central differences
        approximation.  This is less accurate, but requires one fewer
        potential call per iteration.
    gradient : float array
        the true gradient at coords.  If first_order is true and gradient
        is not None then one potential call will be saved. 
    """

    def __init__(self, coords, pot, orthogZeroEigs=0, dx=1e-6,
                 first_order=True, gradient=None, verbosity=1):
        self.pot = pot
        self.first_order = first_order
        self.nfev = 0
        self.verbosity = verbosity
        self.update_coords(coords, gradient=gradient)
        if orthogZeroEigs == 0:
            self.orthogZeroEigs = orthogopt
        else:
            self.orthogZeroEigs = orthogZeroEigs
        self.diff = dx

    def _get_true_energy_gradient(self, coords):
        """return the true energy and gradient at coords"""
        self.nfev += 1
        return self.pot.getEnergyGradient(coords)

    def update_coords(self, coords, gradient=None):
        """update the position at which the curvature is computed"""
        self.coords = coords.copy()
        if self.first_order:
            if gradient is not None:
                # self.true_energy = energy
                self.true_gradient = gradient.copy()
            else:
                if self.verbosity > 1:
                    print("possibly computing gradient unnecessarily")
                true_energy, self.true_gradient = self._get_true_energy_gradient(self.coords)

    def getEnergy(self, vec_in):
        """return the 'energy' a.k.a. the curvature at coords along the direction vec_in"""
        if self.orthogZeroEigs is not None:
            vec_in /= np.linalg.norm(vec_in)
            vec_in = self.orthogZeroEigs(vec_in, self.coords)
        vec = vec_in / np.linalg.norm(vec_in)

        coordsnew = self.coords + self.diff * vec
        Eplus, Gplus = self._get_true_energy_gradient(coordsnew)

        if self.first_order:
            curvature = np.dot((Gplus - self.true_gradient), vec) / self.diff

        else:
            coordsnew = self.coords - self.diff * vec
            Eminus, Gminus = self._get_true_energy_gradient(coordsnew)

            curvature = np.dot((Gplus - Gminus), vec) / (2.0 * self.diff)
        return curvature


    def getEnergyGradient(self, vec_in):
        """return the curvature and the gradient of the curvature w.r.t. vec_in  
        
        vec_in : array 
            A guess for the lowest eigenvector.  It should be normalized
        """
        vecl = 1.
        if self.orthogZeroEigs is not None:
            vec_in /= np.linalg.norm(vec_in)
            vec_in = self.orthogZeroEigs(vec_in, self.coords)
        vec = vec_in / np.linalg.norm(vec_in)

        coordsnew = self.coords + self.diff * vec
        Eplus, Gplus = self._get_true_energy_gradient(coordsnew)
        if self.first_order:
            curvature = np.dot((Gplus - self.true_gradient), vec) / self.diff
        else:
            coordsnew = self.coords - self.diff * vec
            Eminus, Gminus = self._get_true_energy_gradient(coordsnew)
            # use second order central difference method.
            curvature = np.dot((Gplus - Gminus), vec) / (2.0 * self.diff)

        # higher order central differences would be more accurate but it cannot be differentiated analytically
        # DIAG = (EPLUS + EMINUS - 2. * ENERGY) / (self.diff)
        # DIAG3=2*(DIAG-DIAG2/2)
        # C  Although DIAG3 is a more accurate estimate of the diagonal second derivative, it
        # C  cannot be differentiated analytically.

        # compute the analytical derivative of the curvature with respect to vec        
        # GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG2*LOCALV(J1)/VECL**2
        if self.first_order:
            grad = (Gplus - self.true_gradient) * 2.0 / self.diff - 2. * curvature * vec
        else:
            grad = (Gplus - Gminus) / (self.diff * vecl ** 2) - 2.0 * curvature * vec / vecl ** 2
        if self.orthogZeroEigs is not None:
            grad = self.orthogZeroEigs(grad, self.coords)

        # Project out any component of the gradient along vec (which is a unit vector)
        # This is a big improvement for DFTB.
        # js850> grad should already be perpendicular to vec.  this helps with any numerical errors
        grad -= np.dot(grad, vec) * vec

        return curvature, grad


class FindLowestEigenVector(object):
    """A class to compute the lowest eigenvector of the Hessian using Rayleigh-Ritz minimization
    
    Parameters
    ----------
    coords : float array
        the point in space at which to compute the lowest eigenvector
    pot : Potential object
        the potential energy function
    eigenvec0 : float array
        the initial guess for the lowest eigenvector
    orthogZeroEigs : callable
        The function which makes a vector orthogonal to the known
        eigenvectors with zero eigenvalues.  The default assumes global
        translational and rotational symmetry
    dx : float
        the local curvature is approximated using points separated by dx
    first_order : bool
        use the first order forward finite differences approximation for
        the curvature rather than the second order central differences
        approximation.  This is less accurate, but requires one fewer
        potential call per iteration.
    gradient : float array
        the true gradient at coords.  If first_order is true and gradient
        is not None then one potential call will be saved.
    minimizer_kwargs : kwargs
        these kwargs are passed to the optimizer which finds the direction 
        of least curvature
    """

    def __init__(self, coords, pot, eigenvec0=None, orthogZeroEigs=0, dx=1e-6,
                 first_order=True, gradient=None, **minimizer_kwargs):

        self.minimizer_kwargs = minimizer_kwargs

        if eigenvec0 is None:
            # this random vector should be distributed uniformly on a hypersphere.
            eigenvec0 = rotations.vec_random_ndim(coords.shape)
        eigenvec0 = eigenvec0 / np.linalg.norm(eigenvec0)

        # change some default in the minimizer unless manually set
        if "nsteps" not in minimizer_kwargs:
            minimizer_kwargs["nsteps"] = 500
        if "logger" not in minimizer_kwargs:
            minimizer_kwargs["logger"] = logging.getLogger("pele.connect.findTS.leig_quench")

        self.eigpot = LowestEigPot(coords, pot, orthogZeroEigs=orthogZeroEigs, dx=dx,
                                   gradient=gradient,
                                   first_order=first_order)
        self.minimizer = MYLBFGS(eigenvec0, self.eigpot, rel_energy=True,
                                 **self.minimizer_kwargs)

    def stop_criterion_satisfied(self):
        """test if the stop criterion is satisfied"""
        return self.minimizer.stop_criterion_satisfied()

    def update_coords(self, coords, energy=None, gradient=None):
        """update the position at which to compute the eigenvector"""
        self.eigpot.update_coords(coords, gradient=gradient)
        state = self.minimizer.get_state()
        ret = self.get_result()
        self.minimizer = MYLBFGS(ret.eigenvec, self.eigpot, rel_energy=True,
                                 **self.minimizer_kwargs)
        self.minimizer.set_state(state)

    def one_iteration(self):
        """do one iteration of the minimizer"""
        self.minimizer.one_iteration()

    def run(self, niter=None):
        """do niter iterations, or until the stop criterion is satisfied"""
        if niter is None:
            self.minimizer.run()
            return self.get_result()
        else:
            for i in range(niter):
                if self.minimizer.stop_criterion_satisfied():
                    break
                self.one_iteration()
            return self.get_result()

    def get_result(self):
        """return the results object"""
        res = self.minimizer.get_result()
        res.eigenval = res.energy
        res.eigenvec = res.coords / np.linalg.norm(res.coords)
        delattr(res, "energy")
        delattr(res, "coords")
        # res.minimizer_state = self.minimizer.get_state()
        res.nfev = self.eigpot.nfev
        return res


def findLowestEigenVector(coords, pot, eigenvec0=None, H0=None, orthogZeroEigs=0, dx=1e-3,
                          first_order=True, gradient=None,
                          **minimizer_kwargs):
    """Compute the lowest eigenvector of the Hessian using Rayleigh-Ritz minimization

    ***orthogZeroEigs is system dependent, don't forget to set it***

    Parameters
    ----------
    coords :
        the coordinates at which to find the lowest eigenvector
    pot :
        potential object
    eigenvec0 :
        the initial guess for the lowest eigenvector (will be random if not
        passed)
    H0 : float
        the initial guess for the diagonal component of the inverse Hessian
    orthogZeroEigs : callable
        this function makes a vector orthogonal to the known zero
        eigenvectors

            orthogZeroEigs=0  : default behavior, assume translational and
                                rotational symmetry
            orthogZeroEigs=None : the vector is unchanged
    first_order : bool
        use the first order forward finite differences approximation for
        the curvature rather than the second order central differences
        approximation.  This is less accurate, but requires one fewer
        potential call per iteration.
    gradient : float array
        the true gradient at coords.  If first_order is true and gradient
        is not None then one potential call will be saved.
    minimizer_kwargs : 
        any additional keyword arguments are passed to the minimizer
    
    See Also
    --------
    FindTransitionState : uses this class
    """
    minimizer_kwargs = minimizer_kwargs.copy()
    if "iprint" not in minimizer_kwargs:
        minimizer_kwargs["iprint"] = 400
    if "tol" not in minimizer_kwargs:
        minimizer_kwargs["tol"] = 1e-6

    optimizer = FindLowestEigenVector(coords, pot, eigenvec0=eigenvec0, H0=H0, orthogZeroEigs=orthogZeroEigs,
                                      dx=dx, first_order=first_order, gradient=gradient, **minimizer_kwargs)
    result = optimizer.run()
    return result


def analyticalLowestEigenvalue(coords, pot):
    """return the lowest eigenvalue and eigenvector of the hessian computed directly"""
    from pele.utils.hessian import get_sorted_eig

    """for testing"""
    hess = pot.getHessian(coords)
    vals, vecs = get_sorted_eig(hess)

    return vals[0], vecs[:, 0]


#
#
# only testing function below here
#
#


#
#
# def testpot2():
# from pele.potentials.lj import LJ
# import itertools
# pot = LJ()
# a = 1.12 #2.**(1./6.)
# theta = 20./360*np.pi
# coords = [ 0., 0., 0., \
# -a, 0., 0., \
# a*np.cos(theta), a*np.sin(theta), 0. ]
# c = np.reshape(coords, [3,3])
# for i, j in itertools.combinations(range(3), 2):
# r = np.linalg.norm(c[i,:] - c[j,:])
# print i, j, r
#
# def testpot1():
#    from pele.potentials.lj import LJ
#    import itertools
#    pot = LJ()
#    a = 1.12 #2.**(1./6.)
#    theta = 60./360*np.pi
#    coords = [ 0., 0., 0., \
#              -a, 0., 0., \
#              -a/2, a*np.cos(theta), 0., \
#              -a/2, -a*np.cos(theta), 0.1 \
#              ]
#    natoms = len(coords)/3
#    c = np.reshape(coords, [-1,3])
#    for i, j in itertools.combinations(range(natoms), 2):
#        r = np.linalg.norm(c[i,:] - c[j,:])
#        print i, j, r 
#    
#    e, g = pot.getEnergyGradient(coords)
#    print "initial E", e
#    print "initial G", g, np.linalg.norm(g)
#
#    eigpot = LowestEigPot(coords, pot)
#    vec = np.random.rand(len(coords))
#    e, g = eigpot.getEnergyGradient(vec)
#    print "eigenvalue", e 
#    print "eigenvector", g
#    
#    if True:
#        e, g, hess = pot.getEnergyGradientHessian(coords)
#        print "shape hess", np.shape(hess)
#        print "hessian", hess
#        u, v = np.linalg.eig(hess)
#        print "max imag value", np.max(np.abs(u.imag))
#        print "max imag vector", np.max(np.abs(v.imag))
#        u = u.real
#        v = v.real
#        print "eigenvalues", u
#        for i in range(len(u)):
#            print "eigenvalue", u[i], "eigenvector", v[:,i]
#        #find minimum eigenvalue, vector
#        imin = 0
#        umin = 10.
#        for i in range(len(u)):
#            if np.abs(u[i]) < 1e-10: continue
#            if u[i] < umin:
#                umin = u[i]
#                imin = i
#        print "lowest eigenvalue ", umin, imin
#        print "lowest eigenvector", v[:,imin]
#
#    
#    from pele.optimize import lbfgs_py as quench
#    ret = quench(vec, eigpot.getEnergyGradient, iprint=10, tol = 1e-5, maxstep = 1e-3, \
#                 rel_energy = True)
#    print ret
#    
#    print "lowest eigenvalue "
#    print umin, imin
#    print "lowest eigenvector"
#    print v[:,imin]
#    print "now the estimate"
#    print ret[1]
#    print ret[0]
#
#def testpot3():
#    from transition_state_refinement import guesstsLJ
#    pot, coords, coords1, coords2 = guesstsLJ()
#    coordsinit = np.copy(coords)
#
#    eigpot = LowestEigPot(coords, pot)
#    
#    vec = np.random.rand(len(coords))
#    
#    from pele.optimize import lbfgs_py as quench
#    ret = quench(vec, eigpot.getEnergyGradient, iprint=400, tol = 1e-5, maxstep = 1e-3, \
#                 rel_energy = True)
#
#    eigval = ret[1]
#    eigvec = ret[0]
#    print "eigenvalue ", eigval
#    print "eigenvector", eigvec
#
#    if True:
#        e, g, hess = pot.getEnergyGradientHessian(coords)
#        u, v = np.linalg.eig(hess)
#        u = u.real
#        v = v.real
#        print "eigenvalues", sorted(u)
#        #for i in range(len(u)):
#        #    print "eigenvalue", u[i], "eigenvector", v[:,i]
#        #find minimum eigenvalue, vector
#        imin = 0
#        umin = 10.
#        for i in range(len(u)):
#            if np.abs(u[i]) < 1e-10: continue
#            if u[i] < umin:
#                umin = u[i]
#                imin = i
#        #print "lowest eigenvalue ", umin, imin
#        #print "lowest eigenvector", v[:,imin]
#        
#        
#        
#        trueval, truevec = u[imin], v[:,imin]
#        print "analytical lowest eigenvalue", trueval
#        maxdiff = np.max(np.abs(truevec - eigvec))
#        print "maximum difference between estimated and analytical eigenvectors", maxdiff, \
#            np.linalg.norm(eigvec), np.linalg.norm(truevec), np.dot(truevec, eigvec)
#        if True:
#            print eigvec
#            print truevec
#
#
#
#if __name__ == "__main__":
#    #testpot1()
#    testpot3()
#    
    
    

