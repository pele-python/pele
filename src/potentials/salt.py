from math import *
import numpy as np
import salt_
import potential

class salt(potential.potential):
    """binary lennard jones potential with smooth cutoff"""
    def __init__(self):
        # self.natoms = natoms
        pass
    
    def getEnergy(self, coords):
        E = salt_.energy(coords)
        #grad=np.zeros(coords.shape[0], np.float64)
        #E = salt_.gradient(coords, grad)
        #print E,coords[-6:]
        return E

    def getEnergyGradient(self, coords):
        grad=np.zeros(coords.shape[0], np.float64)
        E = salt_.gradient(coords, grad)
        #print E,coords[-6:]
        #print E
        #print coords
        #print grad
        #grad[-6:]=0.
        return E, grad 

    def toReal(self, coords):
        salt_.toReal(coords);