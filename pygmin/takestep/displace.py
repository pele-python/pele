'''
Created on Jun 6, 2012

@author: vr274
'''

import numpy as np

class RandomDisplacement(object):
    def __init__(self, stepsize=0.3, first=None, last=None):
        self.stepsize = stepsize
        self.first = first
        self.last = last

    def takeStep(self, coords, **kwargs):
        c = coords[self.first:self.last]
        c += self.stepsize*(np.random.random(len(c))-0.5)*2.
        pass
        
    def updateStep(self, accepted, **kwargs):
        pass
            
class UniformDisplacement(RandomDisplacement):        
    def takeStep(self, coords, **kwargs):
        import rotations
        c = coords[self.first:self.last]        
        for x in c.reshape(c.size/3,3):
            x += self.stepsize*rotations.vec_random()

class RotationalDisplacement(RandomDisplacement):
    def orientational_step(self, coords, **kwargs):
        """
        take a random orientational step
        """
        import rotations
        c = coords[self.first:self.last]        
        for x in c.reshape(c.size/3,3):
            rotations.takestep_aa( x, self.stepsize ) 

    