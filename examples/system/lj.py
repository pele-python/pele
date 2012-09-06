import numpy as np
from pygmin.potentials.lj import LJ
from pygmin import system

natoms = 13

class LJSystem(system.ClusterSystem):    
    def create_potential(self):
        return LJ()
        
    def initial_coords(self):
        return np.random.random(3*natoms)
    
system.run_basinhopping(LJSystem())