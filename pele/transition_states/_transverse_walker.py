"""routines for minimizing a function in the space perpendicular to a given vector
"""
from __future__ import print_function
import numpy as np

from pele.potentials import BasePotential
from pele.optimize import LBFGS

class _TransversePotential(BasePotential):
    """This wraps a potential and returns the gradient with the component parallel to a vector removed
    
    Parameters
    ----------
    potential : Potential object
    vector : float array
        the vector along which to remove the gradient
    
        
    """
    def __init__(self, potential, vector):
        self.pot = potential
        self.update_vector(vector)
        self.nfev = 0

    def update_vector(self, vector):
        """update the vector"""
        self.vector = vector.copy()
        # normalize
        self.vector /= np.linalg.norm(vector)

    def projected_energy_gradient(self, true_energy, true_gradient):
        """return the projected energy and gradient without computing the gradient"""
        grad = true_gradient - np.dot(true_gradient, self.vector) * self.vector
        return true_energy, grad

    def _get_true_energy_gradient(self, coords):
        """compute the true energy and gradient at coords"""
        self.nfev += 1
        e, grad = self.pot.getEnergyGradient(coords)
        self._true_gradient = grad.copy()
        self._true_energy = e
        self._true_coords = coords.copy()
        return e, grad

    def get_true_energy_gradient(self, coords):
        """return the true energy and gradient at coords
        
        this should primarily be used to access the energy and gradient that have
        already been computed
        """
        if (coords == self._true_coords).all():
            return self._true_energy, self._true_gradient.copy()
        else:
            print("warning: get_true_gradient should only be used to access precomputed energies and gradients")
#            raise Exception("get_true_gradient should only be used to access precomputed energies and gradients")
            return self._get_true_energy_gradient(coords)
    
    def getEnergyGradient(self, coords):
        """return the energy and the gradient with the component along the self.vector removed.  
        
        For use in energy minimization in the space perpendicular to eigenvec
        """
        e, grad = self._get_true_energy_gradient(coords)
        return self.projected_energy_gradient(e, grad)

class _TransverseWalker(object):
    """It minimizes the energy in the direction perpendicular to a vector
    
    this class manages the minimization _TransversePotential
    
    Parameters
    ----------
    coords : float array
        the starting coordinates
    potential : Potential object
    eigenvec : float array
        energy will be minimized in the direction perpendicular to this vector
    energy, gradient : float and float array
        the energy and gradient at position coords
    minimizer_kwargs : kwargs
        these kwargs are passed to the minimizer
    
    """
    def __init__(self, coords, potential, eigenvec, energy=None, gradient=None, **minimizer_kwargs):
        self.tspot = _TransversePotential(potential, eigenvec)
        if energy is not None and gradient is not None:
            transverse_energy, transverse_gradient = self.tspot.projected_energy_gradient(energy, gradient)
        else:
            transverse_energy, transverse_gradient = None, None

        self.walker = LBFGS(coords, self.tspot,
                            energy=transverse_energy, gradient=transverse_gradient,
                            **minimizer_kwargs)
    
    def update_eigenvec(self, eigenvec, eigenval):
        """update the vecotr"""
        self.tspot.update_vector(eigenvec)
    
    def update_maxstep(self, maxstep):
        """update the maximum step size of the minimizer"""
        self.walker.maxstep = float(maxstep)
        
    def update_coords(self, coords, true_energy, true_gradient):
        """update the position of the optimizer
        
        this must be called after update_eigenvec
        """
        energy, gradient = self.tspot.projected_energy_gradient(true_energy, true_gradient)
        self.walker.update_coords(coords, energy, gradient)
    
    def stop_criterion_satisfied(self):
        """test if the stop criterion is satisfied"""
        return self.walker.stop_criterion_satisfied()

    def get_true_energy_gradient(self, coords):
        """return the true energy and gradient"""
        return self.tspot.get_true_energy_gradient(coords)

#    def get_energy(self):
#        """return the true energy
#        
#        warning it's possible for this to return the wrong energy if the minimizer
#        had an aborted line search on the last iteration.
#        """
#        return self.tspot.true_energy
#
#    def get_gradient(self):
#        """return the true gradient
#        
#        warning it's possible for this to return the wrong energy if the minimizer
#        had an aborted line search on the last iteration.
#        """
#        return self.tspot.true_gradient

    def get_result(self):
        """return the results object"""
        ret = self.walker.get_result()
        ret.nfev = self.tspot.nfev
        return ret
    
    def run(self, niter):
        """do a specified number of iterations, or until the stop criterion is satisfied"""
        for i in range(niter):
            if self.stop_criterion_satisfied():
                break
            self.walker.one_iteration()
        return self.get_result()

