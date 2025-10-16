from pele.potentials import BasePotential
from ._mypotential import mypot_energy_gradient

class MyPotFortran(BasePotential):
    """a Lennard Jones potential with altered exponents
    
    V(r) = r**-24 - r**-12
    """
    def __init__(self, natoms):
        self.natoms = natoms #number of atoms
        self.sig = 1.
        self.eps = 1.
    
    def getEnergyGradient(self, coords):
        e, grad = mypot_energy_gradient(coords, self.eps, self.sig)
        return e, grad
    
    def getEnergy(self, coords):
        e, grad = self.getEnergyGradient(coords)
        return e

def test():
    import numpy as np
    natoms = 13
    x = np.random.rand(3*natoms) * 5
    pot = MyPotFortran(natoms)
    pot.test_potential(x)

if __name__ == "__main__":
    test()
