import numpy as np

from pele.potentials import BasePotential
from pele.optimize import LBFGS

class _TransversePotential(BasePotential):
    """This wraps a potential and returns the gradient with the component parallel to a vector removed
    
    Parameters
    ----------
    
    
        
    """
    def __init__(self, potential, vector):
        """
        """
        self.pot = potential
        self.update_vector(vector)
        self.nfev = 0
    
    def update_vector(self, vector):
        self.vector = vector.copy()
        # normalize
        self.vector /= np.linalg.norm(vector)

    def projected_energy_gradient(self, true_energy, true_gradient):
        grad = true_gradient - np.dot(true_gradient, self.vector) * self.vector
        return true_energy, grad
        
    def getEnergyGradient(self, coords):
        """
        return the energy and the gradient with the component along the
        eigenvec removed.  For use in energy minimization in the space
        perpendicular to eigenvec
        """
        self.nfev += 1
        e, grad = self.pot.getEnergyGradient(coords)
        self.true_gradient = grad.copy()
        self.true_energy = e
        #norm = np.sum(self.eigenvec)
        return self.projected_energy_gradient(e, grad)

class _TransverseWalker(object):
    def __init__(self, coords, potential, eigenvec, energy, gradient, **minimizer_kwargs):
        self.tspot = _TransversePotential(potential, eigenvec)
        transverse_energy, transverse_gradient = self.tspot.projected_energy_gradient(energy, gradient) 

        self.walker = LBFGS(coords, self.tspot,
                            energy=transverse_energy, gradient=transverse_gradient,
                            **minimizer_kwargs)
    
    def update_eigenvec(self, eigenvec, eigenval):
        self.tspot.update_vector(eigenvec)
    
    def update_maxstep(self, maxstep):
        self.walker.maxstep = float(maxstep)
        
    def update_coords(self, coords, true_energy, true_gradient):
        """update the position of the optimizer
        
        this must be called after update_eigenvec
        """
        energy, gradient = self.tspot.projected_energy_gradient(true_energy, true_gradient)
        self.walker.update_coords(coords, energy, gradient)
    
    def stop_criterion_satisfied(self):
        return self.walker.stop_criterion_satisfied()

    def get_energy(self):
        """return the true energy"""
        return self.tspot.true_energy

    def get_gradient(self):
        """return the true gradient"""
        return self.tspot.true_gradient

    def get_result(self):
        ret = self.walker.get_result()
        ret.nfev = self.tspot.nfev
        return ret
    
    
    def run(self, niter):
        for i in xrange(niter):
            if self.stop_criterion_satisfied():
                break
            self.walker.one_iteration()
        return self.get_result()

