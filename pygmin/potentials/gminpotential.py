'''
Created on 3 Jun 2012

@author: ruehle
'''

import potential
import numpy as np

class GMINPotential(potential.potential):
    '''
    classdocs
    '''


    def __init__(self, GMIN):
        '''
        Constructor
        '''
        self.GMIN = GMIN
    
    def getEnergy(self, coords):
        return self.GMIN.getEnergy(coords)
        
    def getEnergyGradient(self, coords):
        grad = np.zeros(3*self.GMIN.getNAtoms()) #coords.shape)
        E = self.GMIN.getEnergyGradient(coords, grad)
        return E,grad[0:coords.size]
    
    def getCoords(self):
        coords = np.zeros(self.GMIN.getDOF())
        self.GMIN.getCoords(coords)
        return coords