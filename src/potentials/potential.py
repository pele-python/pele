'''
Created on 13 Apr 2012

@author: ruehle
'''

import numpy as np

class potential(object):
    '''
    Base class for all potentials
    '''
        
    def getEnergyGradient(self, coords):
        return self.getEnergyGradientNumerical(coords)

    def getEnergyGradientNumerical(self, coords):
        return self.getEnergy(coords), self.NumericalDerivative(coords, 1e-8)
            
    def NumericalDerivative(self, coords, eps):
        g = np.zeros(coords.size)
        x = coords.copy()
        for i in xrange(coords.size):
            x[i] += eps
            g[i] = self.getEnergy(x)
            x[i] -= 2. * eps
            g[i] -= self.getEnergy(x)
            g[i] /= 2. * eps
            x[i] += eps
        return g
    
    def getGradient(self, coords):
        e, g = self.getEnergyGradient(coords)
        return g     
