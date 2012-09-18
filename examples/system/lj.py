import numpy as np
from pygmin.potentials.lj import LJ
from pygmin.application import AppClusterBH

natoms = 13

class LJSystem(AppClusterBH):    
    def create_potential(self):
        return LJ()
        
    def initial_coords(self):
        return np.random.random(3*natoms)
    
LJSystem().execute()