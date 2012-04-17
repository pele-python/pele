# -*- coding: iso-8859-1 -*-
import numpy as np


class takeStep:
    """Class to take a monte carlo step.  The stepsize can constant, or
    accessed through function getStep.  RNG must be a random number generator
    returning floats in [0,1)"""
    def __init__(self, RNG = np.random.rand, getStep=None, stepsize=0.3):
        self.RNG = RNG #random number generator
        if getStep == None: 
            self.getStep = lambda : stepsize
        else:
            self.getStep = getStep
        return

    def takeStep(self, coords):
        for j in xrange(len(coords)):
            rand = 2.*self.RNG()-1.
            #print "rand ", rand
            coords[j] += self.getStep()*rand
