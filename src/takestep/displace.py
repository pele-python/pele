'''
Created on Jun 6, 2012

@author: vr274
'''

import numpy as np

class CartesianDisplacement(object):
    '''
    classdocs
    '''

    def __init__(self, stepsize=0.3, RNG = np.random.rand):
        self.RNG = RNG #random number generator
        self.stepsize = stepsize
    
    def takeStep(self, coords):
        coords += self.stepsize*(self.RNG(len(coords))-0.5)*2.
        
    def updateStep(self, accepted):
        pass
