import numpy as np
import stock.stockaa_ as s

class LJ:
    """binary lennard jones potential with smooth cutoff"""
    def __init__(self, mu=0.3):
        print "using lenard jones cpp implementation"
        self.mu = mu

    def getEnergy(self, coords):
        print "getting energy only"
        grad=np.zeros(coords.shape[0], np.float64)
        E = s.gradient(coords, grad, self.mu)
        return E

    def getEnergyGradient(self, coords):
        grad=np.zeros(coords.shape[0], np.float64)
        E = s.gradient(coords, grad, self.mu)
        return E, grad 
