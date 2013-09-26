import numpy as np
from numpy import exp, sqrt, cos, pi, sin

from pele.potentials import BasePotential
from pele.systems import BaseSystem

class Beale(BasePotential):
    target_E = 0.
    target_coords = np.array([3., 0.5])
    xmin = np.array([-4.5, -4.5])
    xmax = np.array([4.5, 4.5])
    def getEnergy(self, coords):
        x, y = coords
        return (1.5 - x + x*y)**2 + (2.25 - x + x * y**2)**2 + (2.625 - x + x * y**3)**2
        
    def getEnergyGradient(self, coords):
        E = self.getEnergy(coords)
        x, y = coords
        dx = 2. * (1.5 - x + x*y) * (-1. + y) + 2. * (2.25 - x + x * y**2) * (-1. + y**2) + 2. * (2.625 - x + x * y**3) * (-1. + y**3)
        dy = 2. * (1.5 - x + x*y) * (x) +       2. * (2.25 - x + x * y**2) * (2. * y * x) + 2. * (2.625 - x + x * y**3) * (3. * x * y**2)
        return E, np.array([dx, dy])
        
class BealeSystem(BaseSystem):
    def get_potential(self):
        return Beale()
    
    def get_random_configuration(self, eps=1e-3):
        pot = self.get_potential()
        xmin, xmax = pot.xmin, pot.xmax
        x = np.random.uniform(xmin[0] + eps, xmax[0] - eps)
        y = np.random.uniform(xmin[1] + eps, xmax[1] - eps)
        return np.array([x,y])
    

if __name__ == "__main__":
    s = BealeSystem()
    f = s.get_potential()
    f.test_potential(f.target_coords)
    print ""
    f.test_potential(s.get_random_configuration())
    f.test_potential(np.array([1.,1.]))#, print_grads=True)
    
#    from base_function import makeplot2d
#    v = 3.
#    makeplot2d(f, nx=30)

