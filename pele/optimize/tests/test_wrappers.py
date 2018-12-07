from __future__ import print_function
import unittest
import numpy as np

from pele.optimize import _quench
from pele.systems import LJCluster

_x0 = np.array([ 0.32384539,  0.5278966 ,  0.55276594,  0.16522439, -0.51696223,
                0.8617451 , -0.42699832,  0.27924851,  1.34690822, -0.66866717,
                1.06990585,  0.54175599,  0.27789836,  0.20221746, -1.23248991,
               -1.02950813, -1.73053047,  1.46386908,  0.25970512,  1.16981168,
                1.40675158,  1.16324999, -0.20594663,  0.51057205, -0.83096288,
                0.6341283 , -1.27923428,  0.95135107,  0.66864409, -0.32954562,
               -1.30558294,  0.51482533, -0.28214219, -0.17866041, -1.978326  ,
               -0.19502835,  0.73588114,  0.14461464,  1.50069807,  0.36819856,
               -0.24505346, -0.22537656,  0.67060212, -1.22846886,  0.14363615,
                0.12835567, -1.10633029, -0.88200269, -0.57019276,  0.02384344,
                0.30595058,  0.05233391,  1.32190053, -1.25120727, -0.9074215 ,
               -0.73272605,  1.00585751, -0.59035945, -0.26828876, -0.73511906,
               -1.14878187, -1.64598547,  0.35324734, -1.4182875 , -0.58012844,
                0.01562395,  0.99903911, -1.09001572,  1.20891914, -0.87801308,
                1.50688001, -0.5457945 ,  1.29998085,  0.8670815 ,  0.83275257,
               -1.61921131, -0.26082059, -1.03603232, -1.01068053, -1.29622377,
               -0.74968498, -0.090193  , -1.58564678,  0.89183383, -0.38380976,
               -0.97007942,  0.06056825, -1.46490078,  0.25993512,  0.82553069,
               -0.16855842,  0.71355274, -0.38843881])

class TestMinimizers(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        natoms = 31
        self.system = LJCluster(natoms)
        self.pot = self.system.get_potential()
        
        # get a partially minimized structure
        self.x0 = _x0.copy()
        self.E0 = -128.289209867
        
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
        res = _quench.bfgs_scipy(self.x0, self.pot, tol=1e-6)
        print("bfgs res", res)
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
        
