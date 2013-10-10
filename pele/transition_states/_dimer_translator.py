"""
tools to invert the gradient along the a given direction and optimize in that space
"""
import numpy as np

from pele.optimize import MYLBFGS, LBFGS

class _DimerTranslator(object):
    """object to manage the translation of the dimer using a optimization algorithm
    """
    def __init__(self, coords, potential, eigenvec, **minimizer_kwargs):
        self.dimer_potential = _DimerPotential(potential, eigenvec)
        self.minimizer = LBFGS(coords, self.dimer_potential, **minimizer_kwargs)

    def stop_criterion_satisfied(self):
        return self.minimizer.stop_criterion_satisfied()

    def get_energy(self):
        """return the true energy"""
        return self.dimer_potential.true_energy

    def get_gradient(self):
        """return the true gradient"""
        return self.dimer_potential.true_gradient
    
    def update_eigenvec(self, eigenvec, eigenval):
        self.dimer_potential.update_eigenvec(eigenvec)
    
    def update_coords(self, coords, true_energy, true_gradient):
        """update the position of the optimizer
        
        this must be called after update_eigenvec
        """
        energy, gradient = self.dimer_potential.projected_energy_gradient(true_energy, true_gradient)
        self.minimizer.update_coords(coords, energy, gradient)

    def update_maxstep(self, maxstep):
        self.minimizer.maxstep = float(maxstep)

    def run(self, niter):
        for i in xrange(niter):
            if self.stop_criterion_satisfied():
                break
            self.minimizer.one_iteration()
        return self.get_result()
    
    def get_result(self):
        return self.minimizer.get_result()
    
    def projected_energy_gradient(self, energy, gradient):
        return self.dimer_potential.projected_energy_gradient(energy, gradient)


class _DimerPotential(object):
    """Wrapper for a Potential object where the gradient is inverted along the direction of the eigenvector
    
    this is used to optimize towards a saddle point
    """
    def __init__(self, potential, eigenvec0,
                  leig_kwargs=None,
                  ):
        self.potential = potential
        self.update_eigenvec(eigenvec0)
        self.nfev = 0
    
    def projected_energy_gradient(self, energy, gradient):
        """return the energy and the gradient at x with the gradient inverted along the eigenvector"""
        projgrad = gradient - 2. * np.dot(gradient, self.eigenvec) * self.eigenvec
#        print "overlap of g and eigenvec", np.dot(gradient, self.eigenvec)
        return 0., projgrad

#    def getEnergyGradientInverted(self, x):
#        """return the energy and the gradient at x with the gradient inverted along the eigenvector"""
#        g -= 2. * np.dot(g, self.eigenvec) * self.eigenvec
#        return e, g

    def update_eigenvec(self, eigenvec):
        self.eigenvec = eigenvec.copy()
        self.eigenvec /= np.linalg.norm(self.eigenvec)
    
    def getEnergyGradient(self, x):
        """gradient at x with the gradient inverted along the eigenvector
        
        the returned energy is 0 because we are not minimizing in the energy 
        """
        e, g = self.potential.getEnergyGradient(x)
        self.true_energy = e
        self.true_gradient = g.copy()
        self.nfev += 1
        return self.projected_energy_gradient(e, g)
    

