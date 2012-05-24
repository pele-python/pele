# -*- coding: iso-8859-1 -*-
import numpy as np
from take_step.adaptive_step import manageStepSize


class takeStep(object):
    """Class to take a monte carlo step.  The stepsize can constant, or
    accessed through function getStep.  RNG must be a random number generator
    returning floats in [0,1)"""
    def __init__(self, stepsize=0.3, RNG = np.random.rand):
        self.RNG = RNG #random number generator
        
        """
        define functions functions:
        self.getStep()
        self.updateStep(accepted)
        """
        self.getStepSize = lambda : stepsize
        self.updateStep = lambda x: x #do nothing

    def useAdaptiveStep(self, stepsize=None, acc_ratio=0.5, freq=100, adaptive_class=None ):
        """
        use adaptive step size
        """
        """
        define functions:
        self.getStep()
        self.updateStep(accepted)
        """
        if stepsize == None:  #use current stepsize
            stepsize = self.getStepSize()
        if adaptive_class != None:
            self.adaptive_class = adaptive_class
            self.getStepSize = self.adaptive_class.getStepSize
            self.updateStep = self.adaptive_class.getStepSize
        else:
            self.adaptive_class = manageStepSize(stepsize, acc_ratio, freq)
            self.getStepSize = self.adaptive_class.getStepSize
            self.updateStep = self.adaptive_class.insertStep
            
    def useFixedStep(self, stepsize=None):
        """
        don't use adaptive step size.  This function can be use to stop
        adaptive step size after a certain number of steps.  E.g. after equilibration
        """
        if stepsize == None:  #use current stepsize
            stepsize = self.getStepSize()
        self.getStepSize = lambda : stepsize
        self.updateStep = lambda x: x #do nothing


            
    def takeStep(self, coords):
        coords += self.getStepSize()*self.RNG(len(coords))            
