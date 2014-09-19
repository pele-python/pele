import unittest
import numpy as np

from pele.optimize import _quench
from pele.systems import LJCluster


class TestMinimizers(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        natoms = 31
        self.system = LJCluster(natoms)
        self.pot = self.system.get_potential()
        
        # get a partially minimized structure
        x0 = self.system.get_random_configuration()
        ret = _quench.lbfgs_py(x0, self.pot, tol=1e-1)
        self.x0 = ret.coords.copy()
        self.E0 = ret.energy
        
        ret = _quench.lbfgs_py(self.x0, self.pot, tol=1e-7)
        self.x = ret.coords.copy()
        self.E = ret.energy
    
    def check_attributes(self, res):
        self.assertTrue(hasattr(res, "energy"))
        self.assertTrue(hasattr(res, "coords"))
        self.assertTrue(hasattr(res, "nsteps"))
        self.assertTrue(hasattr(res, "nfev"))
        self.assertTrue(hasattr(res, "rms"))
        self.assertTrue(hasattr(res, "grad"))
        self.assertTrue(hasattr(res, "success"))
    
    def test_lbfgs_py(self):
        res = _quench.lbfgs_py(self.x0, self.pot, tol=1e-7, debug=True)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
        
    def test_lbfgs_cpp(self):
        res = _quench.lbfgs_cpp(self.x0, self.pot, tol=1e-7)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
        
    def test_mylbfgs(self):
        res = _quench.mylbfgs(self.x0, self.pot, tol=1e-7)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
    
    def test_fire(self):
        res = _quench.fire(self.x0, self.pot, tol=1e-7)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
    
    def test_modifed_fire(self):
        res = _quench.modifiedfire_cpp(self.x0, self.pot, tol=1e-7)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
    
    def test_lbfgs_scipy(self):
        res = _quench.lbfgs_scipy(self.x0, self.pot, tol=1e-7)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
    
    def test_bfgs_scipy(self):
        res = _quench.bfgs_scipy(self.x0, self.pot, tol=1e-7)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 4)
        self.check_attributes(res)
        
    def test_steepest_descent(self):
        res = _quench.steepest_descent(self.x0, self.pot, tol=1e-1)
        self.assertTrue(res.success)
        self.assertAlmostEqual(self.E, res.energy, 1)
        self.check_attributes(res)
        
        
if __name__ == "__main__":
    unittest.main()
        