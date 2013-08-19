import numpy as np

from pele.potentials import BasePotential
from pele.potentials.fortran.morse import morse as fmorse


class Morse(BasePotential):
    """
    the morse potential
    
    R0 is the position of the bottom of the well.
    rho is the width of the well and has units of inverse length.
    A is the energy scale.

    """
    def __init__(self, rho=1.6047, r0=2.8970, A=0.7102):
        self.rho = rho
        self.r0 = r0
        self.A = A
    
    def getEnergy(self, x):
        v, emorse = fmorse(x, False, self.rho, self.r0, self.A)
        return emorse

    def getEnergyGradient(self, x):
        v, emorse = fmorse(x, True, self.rho, self.r0, self.A)
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