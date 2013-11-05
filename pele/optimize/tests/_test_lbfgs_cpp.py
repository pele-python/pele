import unittest
import numpy as np
import os

from pele.potentials._pythonpotential import CppPotentialWrapper
from pele.potentials import BasePotential, _lj_cpp
from pele.optimize import LBFGS_CPP

ndof = 4
_xrand = np.random.uniform(-1,1,[ndof])
_xmin = np.zeros(ndof)
_emin = 0.

class _E(BasePotential):
    def getEnergy(self, x):
        return np.dot(x, x)

class _EG(object):
    def getEnergy(self, x):
        return np.dot(x, x)
    def getEnergyGradient(self, x):
        return self.getEnergy(x), 2. * x

class _Raise(BasePotential):
    def getEnergy(self, x):
        raise NotImplementedError
    def getEnergyGradient(self, x):
        raise NotImplementedError

class TestLBFGS_CPP_PP(unittest.TestCase):
    def test_raises(self):
        with self.assertRaises(NotImplementedError):
            lbfgs = LBFGS_CPP(_xrand, _Raise())
        with self.assertRaises(NotImplementedError):
            lbfgs = LBFGS_CPP(_xrand, _Raise())
            lbfgs.run()

class TestLBFGS_CPP(unittest.TestCase):
    def do_test(self, pot):
        lbfgs = LBFGS_CPP(_xrand, pot)
        res = lbfgs.run()
        self.assertAlmostEqual(res.energy, _emin, 4)
        self.assertTrue(res.success)
        self.assertLess(np.max(np.abs(res.coords - _xmin)), 1e-2)
        self.assertGreater(res.nfev, 0)
        
    def test_E(self):
        self.do_test(_E())

    def test_EG(self):
        self.do_test(_EG())

class TestLBFGS_CPP_Raises(unittest.TestCase):
    def test_raises(self):
        pot = _lj_cpp._ErrorPotential()
        with self.assertRaises(RuntimeError):
            lbfgs = LBFGS_CPP(_xrand, pot)
        with self.assertRaises(RuntimeError):
            lbfgs = LBFGS_CPP(_xrand, pot)
            lbfgs.run()


if __name__ == "__main__":
    unittest.main()
