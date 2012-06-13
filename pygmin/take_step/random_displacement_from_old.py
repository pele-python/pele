# -*- coding: iso-8859-1 -*-
import numpy as np


class takeStep:
    """
    Class to take a monte carlo step from previous, unquenched, coordinates.  
    
    The stepsize can be constant, or accessed through function getStep. 
    
    This take_step routine preserves the canonical ensemble, in contrast to normal basin hopping.  
    Statistics will be thermodynamically meaningful.
    
    note: this step is not rotationally invariant 
    """
    def __init__(self, getStep=None, stepsize=0.3):
        if getStep == None: 
            self.getStep = lambda : stepsize
        else:
            self.getStep = getStep
        return

    def takeStep(self, coords):
        """
        take a random step from the saved coordinates in self.coords_save
        """
        try:
            self.coords_save  #see if self.coords_save is defined
        except:
            self.coords_save = coords.copy() #first step, must initialize coords_save
        coords[:] = self.coords_save[:] + self.getStep()*np.random.uniform(-1,1,len(coords))
        self.coords_new = coords.copy()
            
    
    def checkAccepted(self, quenchedE, quenched_coords, accepted):
        """
        if the step was accepted, update the coordinates in self.coords_save
        """
        if accepted:
            self.coords_save = self.coords_new
        
