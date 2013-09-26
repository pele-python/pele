import numpy as np
import unittest

class _BaseTest(unittest.TestCase):
    
#    def energy_test(self, x, e):
#        e = self.pot.getEnergy(x)
#        self.assertAlmostEqual(e, self.target_E, 4)
    
    def grad_test(self, x):
        e, g = self.pot.getEnergyGradient(x)
        e1, numerical_g = self.pot.getEnergyGradientNumerical(x)
        self.assertLess(np.max(np.abs(g - numerical_g)), 1e-3)
        self.assertAlmostEqual(e, e1, 4)
    
    def test_e_min(self):
        e = self.pot.getEnergy(self.xmin)
        self.assertAlmostEqual(e, self.Emin, 4)
        
    def test_grad_min(self):
        e, g = self.pot.getEnergyGradient(self.xmin)
        self.assertAlmostEqual(e, self.Emin, 4)
        self.assertLess(np.max(np.abs(g)), 1e-4)
        self.grad_test(self.xmin)

    def test_hess_min(self):
        h = self.pot.getHessian(self.xmin)
        eigenvals = np.linalg.eigvals(h)
        self.assertGreater(np.min(eigenvals), 0.)
    
    def test_random(self):
        self.grad_test(self.xrandom)
    

    