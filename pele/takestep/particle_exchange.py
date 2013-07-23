import numpy as np
import random

from pygmin.takestep.generic import Takestep

__all__ = ["ParticleExchange"]

class ParticleExchange(Takestep):
    """Implement a takestep move which swaps two un-like atoms
    
    Choose a random atom from group A and a random atom from group B
    and exchange their xyz coordinates.
    
    Parameters
    ----------
    Alist, Blist : list of integers
        the indices of the atoms in each of the groups
    verbose : bool
        print debugging info
    """
    def __init__(self, Alist, Blist, verbose=False):
        self.Alist = np.array(Alist)
        self.Blist = np.array(Blist)
        self.verbose = verbose
        
        self.naccept = 0
        self.ntry = 0

    def takeStep(self, coords, **kwargs):
        iA = random.choice(self.Alist)
        iB = random.choice(self.Blist)
        if self.verbose:
            print "exchange atoms", iA, iB, "accepted", self.naccept, "out of", self.ntry
        
        coords = coords.reshape(-1,3)
        temp = coords[iA,:].copy()
        coords[iA,:] = coords[iB,:]
        coords[iB,:] = temp
        self.ntry += 1
        return coords
    
    def updateStep(self, accepted, **kwargs):
        '''feedback from basin hopping if last step was accepted'''
        if accepted:
            self.naccept += 1
        


import unittest

class TestParticleExchange(unittest.TestCase):
    def setUp(self):
        Alist = [1,3,5]
        Blist = [0,2,4]
        x = np.zeros([6,3], int)
        for i in Alist:
            x[i,:] = [1,1,1]
        
        self.step = ParticleExchange(Alist, Blist)
        self.x = x.flatten()
    
    def test(self):
        x = self.x
        self.step(x)
        x = x.reshape(-1,3)
        sA = x[self.step.Alist,:].sum()
        self.assertEqual(sA, 6)
        sB = x[self.step.Blist,:].sum()
        self.assertEqual(sB, 3)

if __name__ == "__main__":
    unittest.main()
        