'''
Created on Jun 6, 2012

@author: vr274
'''

from .generic import TakestepInterface
import numpy as np

class AdaptiveStepsize(TakestepInterface):
    '''
    Adaptive stepsize
    '''

    def __init__(self, stepclass, acc_ratio=0.5, factor=0.9, frequency=1000):
        self.stepclass = stepclass
        self.accrat = acc_ratio #target accept ratio            
        self.factor = factor
        self.nstepsaccrat = frequency

        self.naccepted = 0
        self.nsteps = 0

    
    def takeStep(self, coords, **kwargs):
        self.stepclass.takeStep(coords)
        
    def updateStep(self, accepted, **kwargs):
        """tell us whether a step was accepted or rejected"""
        self.nsteps += 1
        if accepted: 
            self.naccepted += 1
        if self.nsteps == self.nstepsaccrat:
            self.adjustStep()
            self.nsteps = 0
                        
    def adjustStep(self):
        """adjust stepsize"""
        rat = float(self.naccepted)/self.nsteps
        if rat > self.accrat:
            self.stepclass.scale( 1./ self.factor)
        else:
            self.stepclass.scale( self.factor )

        self.nsteps = 0
        self.naccepted = 0
        print "accrat was ", rat, "new stepsize is ", self.stepclass.stepsize, "f is", self.factor        