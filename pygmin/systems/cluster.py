import numpy as np

from pygmin.systems import BaseSystem, NotImplemented
from pygmin.potentials import LJ
from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
from pygmin.transition_states import orthogopt
from pygmin.mindist import minPermDistStochastic, MinDistWrapper, ExactMatchCluster
from compiler.ast import Not

__all__ = ["AtomicCluster"]

class AtomicCluster(BaseSystem):
    """
    Define an atomic cluster.  
    
    This is a system of point particles with global rotational 
    and translational symmetry and some form of permutational 
    symmetry.
    """
    def get_potential(self):
        return LJ(self.natoms)
    
    def get_random_configuration(self):
        coords = np.random.uniform(-1, 1, [3*self.natoms]) * 1.5 * float(self.natoms)**(1./3)
        return coords
    
    def get_permlist(self):
        raise NotImplemented
    
    def get_compare_exact(self, **kwargs):
        """this function quickly determines whether two clusters are identical
        given translational, rotational and permutational symmeties
        """
        raise NotImplemented
        permlist = self.get_permlist()
        return ExactMatchCluster(permlist=permlist, **kwargs)
    
    def get_mindist(self, **kwargs):
        """return a function which puts two structures in best alignment.
        
        take into account global rotational symmetry, global translational
        symmetry and permutational symmetry
        """
        permlist = self.get_permlist()
        return MinDistWrapper(minPermDistStochastic, permlist=permlist, **kwargs)
        
    def get_orthogonalize_to_zero_eigenvectors(self):
        """the zero eigenvectors correspond to 3 global translational
        degrees of freedom and 3 global rotational degrees of freedom"""
        return orthogopt

    def get_takestep(self, stepsize=0.6, **kwargs):
        """random displacement with adaptive step size 
        adaptive temperature"""
        takeStep = RandomDisplacement(stepsize=stepsize)
        tsAdaptive = AdaptiveStepsizeTemperature(takeStep, **kwargs)
        return tsAdaptive

