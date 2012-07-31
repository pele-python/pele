'''
Created on Jun 6, 2012

@author: vr274
'''

import numpy as np
from generic import TakestepSlice
from pygmin.utils import rotations

class RandomDisplacement(TakestepSlice):
    def __init__(self, stepsize=1.0):
        TakestepSlice.__init__(self, stepsize=stepsize)
    def takeStep(self, coords, **kwargs):
        c = coords[self.srange]
        c += self.stepsize*(np.random.random(c.shape)-0.5)*2.
            
class UniformDisplacement(TakestepSlice):        
    def takeStep(self, coords, **kwargs):
        c = coords[self.srange]        
        for x in c.reshape(c.size/3,3):
            x += self.stepsize*rotations.vec_random()

class RotationalDisplacement(TakestepSlice):
    def takeStep(self, coords, **kwargs):
        """
        take a random orientational step
        """
        c = coords[slice]        
        for x in c.reshape(c.size/3,3):
            rotations.takestep_aa( x, self.stepsize ) 

    