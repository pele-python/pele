import numpy as np

from pele.potentials import BasePotential
from pele.potentials.fortran.morse import morse as fmorse


class Morse(BasePotential):
    def __init__(self, rho=5.):
        self.rho = rho
    
    def getEnergy(self, x):
        v, emorse = fmorse(x, False, self.rho)
        return emorse

    def getEnergyGradient(self, x):
        v, emorse = fmorse(x, True, self.rho)
        return emorse, v


def test():
    m = Morse(1.)
    x = np.random.rand(3*10)
    e = m.getEnergy(x)
    print e
    
    e, v = m.getEnergyGradient(x)
    print e, v

if __name__ == "__main__":
    test()