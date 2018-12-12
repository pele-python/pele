"""
tools to invert the gradient along the a given direction and optimize in that space
"""
from __future__ import print_function
import numpy as np

from pele.optimize import LBFGS

class _DimerTranslator(object):
    """object to manage the translation of the dimer using an optimization algorithm
    
    Parameters
    ----------
    coords : float array
        the starting point of the dimer
    potential : Potential object
    eigenvec : float array
        the initial direction along which the dimer lies
    minimizer_kwargs : kwargs
        these kwargs are passed to the optimizer
    """
    def __init__(self, coords, potential, eigenvec, **minimizer_kwargs):
        self.dimer_potential = _DimerPotential(potential, eigenvec)
        self.minimizer = LBFGS(coords, self.dimer_potential, **minimizer_kwargs)

    def stop_criterion_satisfied(self):
        """test if the stop criterion is satisfied"""
        return self.minimizer.stop_criterion_satisfied()

    def get_true_energy_gradient(self, coords):
        """return the true energy and gradient"""
        return self.dimer_potential.get_true_energy_gradient(coords)

#    def get_energy(self):
#        """return the true energy"""
#        return self.dimer_potential.true_energy
#
#    def get_gradient(self):
#        """return the true gradient"""
#        return self.dimer_potential.true_gradient
    
    def update_eigenvec(self, eigenvec, eigenval):
        """update the direction (rotation) of the dimer"""
        self.dimer_potential.update_eigenvec(eigenvec)
    
    def update_coords(self, coords, true_energy, true_gradient):
        """update the position of the dimer
        
        this must be called after update_eigenvec
        """
        energy, gradient = self.dimer_potential.projected_energy_gradient(true_energy, true_gradient)
        self.minimizer.update_coords(coords, energy, gradient)

    def update_maxstep(self, maxstep):
        """change the maximum step size of the optimizer"""
        self.minimizer.maxstep = float(maxstep)

    def run(self, niter):
        """do a specified number of iterations, or until the stop criterion is satisfied"""
        for i in range(niter):
            if self.stop_criterion_satisfied():
                break
            self.minimizer.one_iteration()
        return self.get_result()
    
    def get_result(self):
        """return the results object"""
        return self.minimizer.get_result()
    
    def projected_energy_gradient(self, energy, gradient):
        """return the projected energy and gradient"""
        return self.dimer_potential.projected_energy_gradient(energy, gradient)


class _DimerPotential(object):
    """Wrapper for a Potential object where the gradient is inverted along the direction of the eigenvector
    
    this is used to optimize towards a saddle point
    """
    def __init__(self, potential, eigenvec0,
                 leig_kwargs=None):
        self.potential = potential
        self.update_eigenvec(eigenvec0)
        self.nfev = 0
    
    def projected_energy_gradient(self, energy, gradient):
        """return the energy and the gradient with the gradient inverted along the eigenvector"""
        projgrad = gradient - 2. * np.dot(gradient, self.eigenvec) * self.eigenvec
#        print "overlap of g and eigenvec", np.dot(gradient, self.eigenvec)
        return 0., projgrad

    def update_eigenvec(self, eigenvec):
        """update the direction (rotation) of the dimer"""
        self.eigenvec = eigenvec.copy()
        self.eigenvec /= np.linalg.norm(self.eigenvec)
    
    def _get_true_energy_gradient(self, coords):
        """compute the true energy and gradient at coords"""
        self.nfev += 1
        e, grad = self.potential.getEnergyGradient(coords)
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

    
    def getEnergyGradient(self, x):
        """return the energy and gradient at x with the gradient along the eigenvector inverted
        
        Notes
        -----
        The returned energy is 0 because we are not minimizing in the energy.  The
        true energy and gradient are stored at each call
        """
        e, g = self._get_true_energy_gradient(x)
        return self.projected_energy_gradient(e, g)
    

