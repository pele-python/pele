import numpy as np

from pygmin.systems import BaseSystem
from pygmin.potentials import LJ
from pygmin.transition_states import orthogopt
from pygmin.mindist import minPermDistStochastic, MinDistWrapper, ExactMatchCluster
from pygmin.landscape import smoothPath
from pygmin.transition_states import NEB, InterpolatedPathDensity


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
        coords = np.random.uniform(-1, 1, [3*self.natoms]) * 1.1 * float(self.natoms)**(1./3)
        return coords
    
    def get_permlist(self):
        raise NotImplementedError
    
    def get_compare_exact(self, **kwargs):
        """this function quickly determines whether two clusters are identical
        given translational, rotational and permutational symmeties
        """
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

    #
    #below here only stuff for the gui
    #

    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smoothPath(path, mindist, **kwargs)
        
    def createNEB(self, coords1, coords2):
        pot = self.get_potential()
        dist = np.linalg.norm(coords1- coords2)
        if dist < 1.: dist = 1
        image_density = 15.
        
        path = InterpolatedPathDensity(coords1, coords2, 
                                       distance=dist, density=image_density)
        return NEB(path, pot)

