import unittest
import numpy as np
import os

import _base_test
from pele.potentials._pythonpotential import CppPotentialWrapper
from pele.potentials import BasePotential

ndof = 4
_xrand = np.random.uniform(-1,1,[ndof])
_xmin = np.zeros(ndof)
_emin = 0.

class _Eonly(BasePotential):
    def getEnergy(self, x):
#         print "getting energy"
        return np.dot(x, x)

class _EG(_Eonly):
    def getEnergyGradient(self, x):
        print "in python getEnergyGradient"
        return self.getEnergy(x), 2. * x

class _Raise(BasePotential):
    def getEnergy(self, x):
        raise NotImplementedError


class TestEonly(_base_test._BaseTest):
    def setUp(self, ndof = 4):
        self.pot = _Eonly()
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin
    
class TestEonlyWrapped(_base_test._BaseTest):
    def setUp(self, ndof = 4):
        self.pot = CppPotentialWrapper(_Eonly())
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin
    
class TestEG(_base_test._BaseTest):
    def setUp(self, ndof = 4):
        self.pot = _EG()
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin
    
class TestEGWrapped(_base_test._BaseTest):
    def setUp(self, ndof = 4):
        self.pot = CppPotentialWrapper(_EG())
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin

class TestRaised(unittest.TestCase):
    """test the raised errors are transmitted properly through c++"""
    def setUp(self):
        self.xrandom = _xrand
    
    def test_raises(self):
        p = _Raise()
        with self.assertRaises(NotImplementedError):
            p.getEnergy(_xrand)
        with self.assertRaises(NotImplementedError):
            p.getEnergyGradient(_xrand)

    def test_raises_wrap(self):
        p = CppPotentialWrapper(_Raise())
        with self.assertRaises(NotImplementedError):
            p.getEnergy(_xrand)
        with self.assertRaises(NotImplementedError):
            p.getEnergyGradient(_xrand)


def simplertest():
    pot = CppPotentialWrapper(_Eonly())
    e = pot.getEnergy(_xrand)
    print "energy", e
    e, g = pot.getEnergyGradient(_xrand)
    print "energy", e
    print "grad", g
    print "hess", pot.NumericalHessian(_xrand)

    print pot.NumericalDerivative(_xrand)
    print "done done done"

if __name__ == "__main__":
#     simplertest()
    unittest.main()
