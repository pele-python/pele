import unittest
import numpy as np

from pele.optimize import _quench_obsolete
from pele.systems import LJCluster

class TestObsoleteMinimizers(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        natoms = 31
        self.system = LJCluster(natoms)
        self.pot = self.system.get_potential()
        
        # get a partially minimized structure
        x0 = self.system.get_random_configuration()
        ret = _quench_obsolete.lbfgs_py(x0, self.pot.getEnergyGradient, tol=1.e-1)
        self.x0 = ret[0]
        self.E0 = ret[1]
        
        ret = _quench_obsolete.lbfgs_py(self.x0, self.pot.getEnergyGradient, tol=1e-7)
        self.x = ret[0]
        self.E = ret[1]
    
    def test_lbfgs_py(self):
        res = _quench_obsolete.lbfgs_py(self.x0, self.pot.getEnergyGradient)
        self.assertAlmostEqual(self.E, res[1], 4)
        
    def test_mylbfgs(self):
        res = _quench_obsolete.mylbfgs(self.x0, self.pot.getEnergyGradient)
        self.assertAlmostEqual(self.E, res[1], 4)
    
    def test_fire(self):
        res = _quench_obsolete.fire(self.x0, self.pot.getEnergyGradient, tol=1e-7)
        self.assertAlmostEqual(self.E, res[1], 4)
    
    def test_lbfgs_scipy(self):
        res = _quench_obsolete.lbfgs_scipy(self.x0, self.pot.getEnergyGradient, tol=1e-7)
        self.assertAlmostEqual(self.E, res[1], 4)
    
    def test_bfgs_scipy(self):
        res = _quench_obsolete.bfgs(self.x0, self.pot.getEnergyGradient, tol=1e-7)
        self.assertAlmostEqual(self.E, res[1], 4)
        
        
if __name__ == "__main__":
    unittest.main()

