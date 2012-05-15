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


    #routines involving interaction lists
    def getEnergyListSlow(self, coords, ilist):
        coords2 = np.zeros(6)
        E = 0.
        for i,j  in ilist:
            coords2[:3] = coords[i*3:i*3+3]
            coords2[3:] = coords[j*3:j*3+3]
            E += self.getEnergy( coords2 )
        return E

    def getEnergyList(self, coords, ilist):
        """
        the energy of select interactions defined in ilist
        """
        return self.getEnergyListSlow(coords, ilist)

    def getEnergyGradientListSlow(self, coords, ilist):
        coords2 = np.zeros(6)
        E = 0.
        g = np.zeros( len(coords) )
        for i,j  in ilist:
            coords2[:3] = coords[i*3:i*3+3]
            coords2[3:] = coords[j*3:j*3+3]
            de, dg = self.getEnergyGradient( coords2 )
            E += de
            g[i*3:i*3+3] += dg[:3]
            g[j*3:j*3+3] += dg[3:]
        return E, g

    def getEnergyGradientList(self, coords, ilist):
        """
        the energy and gradient of select interactions defined in ilist
        """
        return self.getEnergyGradientListSlow(coords, ilist)
