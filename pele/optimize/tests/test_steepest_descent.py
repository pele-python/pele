# see: test_lbfgs_cpp.py

from __future__ import division

import numpy as np
import os
import unittest as ut

from pele.optimize import SteepestDescentCPP
from pele.potentials import BasePotential
from pele.potentials import LJ

ndof = 4
_xrand = np.random.uniform(-1, 1, [ndof])
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
    
class TestSDCPP(ut.TestCase):
    def test_raises(self):
        with self.assertRaises(NotImplementedError):
            gd = SteepestDescentCPP(_xrand, _Raise())
            gd.run()
            
class TestSDCPP(ut.TestCase):
    def do_check(self, pot):
        gd = SteepestDescentCPP(_xrand, pot)
        res = gd.run()
        self.assertAlmostEqual(res.energy, _emin, 4)
        self.assertTrue(res.success)
        self.assertLess(np.max(np.abs(res.coords - _xmin)), 1e-2)
        self.assertGreater(res.nfev, 0)

    def test_E(self):
        self.do_check(_E())

    def test_EG(self):
        self.do_check(_EG())

    def assert_same(self, res1, res2):
        self.assertEqual(res1.energy, res2.energy)
        self.assertEqual(res1.rms, res2.rms)
        self.assertEqual(res1.nfev, res2.nfev)

    def test_run_niter(self):
        gd1 = SteepestDescentCPP(_xrand, _EG())
        res1 = gd1.run()
        gd2 = SteepestDescentCPP(_xrand, _EG())
        res2 = gd2.run(res1.nsteps)
        self.assert_same(res1, res2)

    def test_run_niter2(self):
        gd1 = SteepestDescentCPP(_xrand, _EG())
        res1 = gd1.run()
        gd2 = SteepestDescentCPP(_xrand, _EG())
        res2 = gd2.run(res1.nsteps / 2)
        res2 = gd2.run()
        self.assert_same(res1, res2)

    def test_run_niter3(self):
        gd1 = SteepestDescentCPP(_xrand, _EG())
        res1 = gd1.run(10)
        gd2 = SteepestDescentCPP(_xrand, _EG())
        res2 = gd2.run(5)
        res2 = gd2.run(5)
        self.assert_same(res1, res2)

if __name__ == "__main__":
    ut.main()
