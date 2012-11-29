import numpy as np
import aautils

class RigidFragment(aautils.AASiteType):
    ''' defines a single rigid fragment 
    
    In the most simple case, this is just a whole molecule
    '''
    
    def __init__(self):
        aautils.AASiteType.__init__(self)
        self.atom_positions = []
        self.atom_types = []
        self.atom_masses = []
        
    def add_atom(self, atomtype, pos, mass=1.0):
        '''Add a new atom to the rigid fragment
        
        Parameters
        ----------
        type: string
            type identifier
        pos: np.array
            position of the atom
        mass: mass of the atom
        '''
        self.atom_types.append(atomtype)
        self.atom_positions.append(pos)
        self.atom_masses.append(mass)
        
    def finalize_setup(self):
        '''finalize setup after all sites have been added
        
        This will shift the center of mass to the origin and calculate
        the total mass and weighted tensor of gyration
        '''
        
        # first correct for the center of mass
        com = np.average(self.atom_positions, axis=0, weights=self.atom_masses)
        for x in self.atom_positions:
            x[:] -= com
        
        # calculate total mass
        self.M = np.sum(self.atom_masses)
        
        # now calculate the weighted moment of inertia tensor
        self.S[:] = 0.
        for x, m in zip(self.atom_positions, self.atom_masses):
            self.S[:] += m*np.outer(x, x)
            
class RBSystem(aautils.AASystem):
    pass