from __future__ import print_function
from __future__ import absolute_import
import unittest

import numpy as np
from pele.potentials._pythonpotential import CppPotentialWrapper, _TestingCppPotentialWrapper

from . import _base_test
from pele.potentials import BasePotential


ndof = 4
_xrand = np.random.uniform(-1, 1, [ndof])
_xmin = np.zeros(ndof)
_emin = 0.


class _Eonly(BasePotential):
    def getEnergy(self, x):
        # print "getting energy"
        return np.dot(x, x)


class _EG(_Eonly):
    def getEnergyGradient(self, x):
        # print "in python getEnergyGradient"
        return self.getEnergy(x), 2. * x


class _Raise(BasePotential):
    def getEnergy(self, x):
        raise NotImplementedError


class _ReturnBad1(BasePotential):
    def getEnergy(self, x):
        return "not a double"

    def getEnergyGradient(self, x):
        return 1., "not a np array"


class _ReturnBad2(BasePotential):
    def getEnergy(self, x):
        return 1.

    def getEnergyGradient(self, x):
        return 1.


class TestEonly(_base_test._BaseTest):
    def setUp(self, ndof=4):
        self.pot = _Eonly()
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin


class TestEonlyWrapped(_base_test._BaseTest):
    def setUp(self, ndof=4):
        self.pot = CppPotentialWrapper(_Eonly())
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin


class TestEG(_base_test._BaseTest):
    def setUp(self, ndof=4):
        self.pot = _EG()
        self.xrandom = _xrand
        self.xmin = _xmin
        self.Emin = _emin


class TestEGWrapped(_base_test._BaseTest):
    def setUp(self, ndof=4):
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


class TestBadReturn(unittest.TestCase):
    """
    test getEnergyGradient returns invalid energy
    test getEnergyGradient returns 
    test getEnergyGradient returns 
    test getEnergyGradient returns 
    """

    def test1(self):
        """test getEnergy returns not double"""
        p = _TestingCppPotentialWrapper(_ReturnBad1())
        with self.assertRaises(TypeError):
            e = p.cpp_get_energy(_xrand)

    def test2(self):
        """test getEnergyGradient returns not tuple"""
        p1 = _ReturnBad1()
        p1.getEnergyGradient = lambda s, x: None
        p = _TestingCppPotentialWrapper(p1)
        with self.assertRaises(TypeError):
            e = p.cpp_get_energy_gradient(_xrand)

    def test3(self):
        """test getEnergyGradient returns invalid energy"""
        p1 = _ReturnBad1()
        p1.getEnergyGradient = lambda x: ("not energy", np.zeros(x.size))
        p = _TestingCppPotentialWrapper(p1)
        with self.assertRaises(TypeError):
            e = p.cpp_get_energy_gradient(_xrand)

    def test4(self):
        """test getEnergyGradient returns invalid grad"""
        p1 = _ReturnBad1()
        p1.getEnergyGradient = lambda x: (1., (1., "not numpy array"))
        p = _TestingCppPotentialWrapper(p1)
        with self.assertRaises(ValueError):
            e = p.cpp_get_energy_gradient(_xrand)

    def test5(self):
        """test getEnergyGradient returns list for grad"""
        p1 = _ReturnBad1()
        p1.getEnergyGradient = lambda x: (1., list(np.zeros(x.size)))
        p = _TestingCppPotentialWrapper(p1)
        # this should not raise anything
        e = p.cpp_get_energy_gradient(_xrand)

    def test6(self):
        """test getEnergyGradient returns grad of wrong size"""
        p1 = _ReturnBad1()
        p1.getEnergyGradient = lambda x: (1., np.zeros(x.size - 1))
        p = _TestingCppPotentialWrapper(p1)
        with self.assertRaises(IndexError):
            e = p.cpp_get_energy_gradient(_xrand)


def simplertest():
    pot = CppPotentialWrapper(_Eonly())
    e = pot.getEnergy(_xrand)
    print("energy", e)
    e, g = pot.getEnergyGradient(_xrand)
    print("energy", e)
    print("grad", g)
    print("hess", pot.NumericalHessian(_xrand))

    print(pot.NumericalDerivative(_xrand))
    print("done done done")


def exceptiontest():
    pot = CppPotentialWrapper(_Raise())
    # e = pot.getEnergy(_xrand)
    # print "energy", e
    # e, g = pot.getEnergyGradient(_xrand)
    # print "energy", e
    #     print "grad", g
    print(pot.NumericalDerivative(_xrand))
    print("hess", pot.NumericalHessian(_xrand))

    print("done done done")


if __name__ == "__main__":
    # simplertest()
    # exceptiontest()
    unittest.main()

