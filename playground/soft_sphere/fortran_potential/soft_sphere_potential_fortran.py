import numpy as np

from pele.potentials import BasePotential
import _soft_sphere

class MyPotFortran(BasePotential):
    """a Lennard Jones potential with altered exponents
    
    V(r) = r**-24 - r**-12
    """
    def __init__(self, boxvec, radii, power=2.5):
        self.boxvec = np.array(boxvec)
        self.radii = np.array(radii, dtype=float)
        self.power = float(power)
    
    def _call_fortran(self, x, compute_hess=False):
        grad, hess, e = _soft_sphere.soft_sphere(x, 
                                               self.boxvec[0],
                                               self.boxvec[1],
                                               self.boxvec[2],
                                               self.radii,
                                               self.power,
                                               True,
                                               compute_hess)
        return e, grad, hess
    def getEnergyGradient(self, x):
        e, grad, hess = self._call_fortran(x, compute_hess=False)
        return e, grad
    def getEnergyGradientHessian(self, x):
        e, grad, hess = self._call_fortran(x, compute_hess=True)
        return e, grad, hess
    
    def getEnergy(self, coords):
        e, grad = self.getEnergyGradient(coords)
        return e


from pele.potentials.tests import _base_test

class TestPot(_base_test._TestConfiguration):
    def setUp(self):
        natoms = 20
        radii = np.random.uniform(.9, 1.2, natoms)
        L = 3
        boxvec = np.ones(3) * L
        x = np.random.uniform(0,L, 3*natoms)
        pot = MyPotFortran(boxvec, radii)
#        pot.test_potential(x)
        
        self.x0 = x
        self.pot = pot
        self.e0 = pot.getEnergy(x)

def test():
    natoms = 20
    radii = np.random.uniform(.9, 1.2, natoms)
    L = 3
    boxvec = np.ones(3) * L
    x = np.random.uniform(0,L, 3*natoms)
    pot = MyPotFortran(boxvec, radii)
    pot.test_potential(x)
    
    e, g, hess = pot.getEnergyGradientHessian(x)
    hnum = pot.NumericalHessian(x)
    print hess
    print hnum
    print hess.shape, hnum.shape
    print hess[0,0]
    print hnum[0,0]


if __name__ == "__main__":
#    test()
    import unittest
    unittest.main()