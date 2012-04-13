'''
Created on 13 Apr 2012

@author: ruehle
'''

import potential
import pygmin
import numpy as np

class PatchyParticle(potential.potential):
    '''
    classdocs
    '''


    def __init__(selfparams):
        '''
        Constructor
        '''
    
    def getEnergy(self, coords):
        natoms = coords.size/3 - 3
        m = self.getLatticeMatrix(coords)
        x = np.dot(coords.reshape([natoms,3], m)).reshape(3*natoms + 3)
        x[-3:-1] = coords[-3:-1]
        return pygmin.pap_energy(coords)
    
    def getEnergyGradient(self, coords):
        natoms = coords.size/3 - 3
        m = self.getLatticeMatrix(coords)
        x = np.dot(coords.reshape([natoms,3], m)).reshape(3*natoms + 3)
        x[-3:-1] = coords[-3:-1]           
        grad = np.zeros(x.shape)
        E = pygmin.pap_gradient(x, grad)
        return E, grad
    
    def getLatticeMatrix(self, coords):
        m = np.zeros([3,3])
        m[0][0] = coords[-3]       
        m[1][1] = coords[-2]
        m[2][2] = coords[-1]
        
    def setLatticeMatrix(self, coords, m):
        coords[-3] = m[0][0]
        coords[-2] = m[1][1]
        coords[-1] = m[2][2]