import unittest
import numpy as np

from pele.potentials import test_functions
from pele.optimize import _quench

class TestOptimizersBeale(unittest.TestCase):
    def setUp(self):
#        np.random.seed(0)
        self.system = test_functions.BealeSystem()
        self.pot = self.system.get_potential()
        
        
        self.minimizers = [_quench.mylbfgs, _quench.lbfgs_py,
                           _quench.bfgs_scipy, _quench.fire
                           ]
        self.x = self.pot.target_coords.copy() 
        self.x += np.random.uniform(-0.2, 0.2, self.x.shape)

    def do_test(self, minimizer, **kwargs):
        ret = minimizer(self.x, self.pot, **kwargs)
        self.assertTrue(ret.success)
        self.assertAlmostEqual(ret.energy, self.pot.target_E, 3)
        self.assertLess(np.max(np.abs(ret.coords - self.pot.target_coords)), 1e-3)
    
    def test_lbfgs_py(self):
        self.do_test(_quench.lbfgs_py)

    def test_mylbfgs(self):
        self.do_test(_quench.mylbfgs)

    def test_fire(self):
        self.do_test(_quench.fire, tol=1e-7)

    def test_lbfgs_scipy(self):
        self.do_test(_quench.lbfgs_scipy)

    def test_bfgs_scipy(self):
        self.do_test(_quench.bfgs_scipy)




if __name__ == "__main__":
    unittest.main()