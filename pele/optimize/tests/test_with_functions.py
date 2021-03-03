from __future__ import print_function
import unittest
import numpy as np

from pele.potentials import test_functions
from pele.optimize import _quench


class TestOptimizersBeale(unittest.TestCase):
    def setUp(self):
        # np.random.seed(0)
        self.system = test_functions.BealeSystem()
        self.pot = self.system.get_potential()

        self.x = self.pot.target_coords.copy()
        self.x += np.random.uniform(-0.2, 0.2, self.x.shape)

    def do_check(self, minimizer, **kwargs):
        ret = minimizer(self.x, self.pot, **kwargs)
        self.assertTrue(ret.success)
        self.assertAlmostEqual(ret.energy, self.pot.target_E, 3)
        self.assertLess(np.max(np.abs(ret.coords - self.pot.target_coords)), 1e-3)

    def test_lbfgs_py(self):
        self.do_check(_quench.lbfgs_py)

    def test_lbfgs_cpp(self):
        self.do_check(_quench.lbfgs_cpp)

    def test_mylbfgs(self):
        self.do_check(_quench.mylbfgs)

    def test_fire(self):
        self.do_check(_quench.fire, tol=1e-7)

    def test_lbfgs_scipy(self):
        self.do_check(_quench.lbfgs_scipy)

    def test_bfgs_scipy(self):
        self.do_check(_quench.bfgs_scipy)


class TestOptimizeBooth(TestOptimizersBeale):
    def setUp(self):
        # np.random.seed(0)
        self.system = test_functions.BoothSystem()
        self.pot = self.system.get_potential()

        self.x = self.pot.target_coords.copy()
        # self.x += np.random.uniform(-0.1, 0.1, self.x.shape)
        self.x = self.system.get_random_configuration()


def mytest():
    system = test_functions.BealeSystem()
    print("do pot")
    pot = system.get_potential()
    print("done pot")
    x = pot.target_coords.copy()
    x += np.random.uniform(-0.2, 0.2, x.shape)
    from pele.optimize import LBFGS_CPP

    lbfgs = LBFGS_CPP(x, pot, verbosity=100)
    print("done setting up")
    lbfgs.run()
    res = _quench.lbfgs_cpp(x, pot, verbosity=100)
    # print res
    print(res)


if __name__ == "__main__":
    # mytest()
    # print "totally finished!!!"
    unittest.main()
