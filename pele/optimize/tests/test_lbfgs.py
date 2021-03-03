from __future__ import print_function
import unittest


import numpy as np

from pele.optimize import LBFGS, MYLBFGS
from pele.systems import LJCluster
from pele.potentials import BasePotential


class DiscontinuousHarmonic(BasePotential):
    def getEnergy(self, x):
        e = np.dot(x, x)
        if x[0] < -1:
            e -= 1
        return e

    def getEnergyGradient(self, x):
        e = self.getEnergy(x)
        g = 2. * x
        return e, g


def arrays_nearly_equal(self, a1, a2, **kwargs):
    if len(kwargs) == 0:
        kwargs = dict(places=5)
    self.assertEqual(a1.shape, a2.shape)
    for v1, v2 in zip(a1.reshape(-1), a2.reshape(-1)):
        self.assertAlmostEqual(v1, v2, **kwargs)


class TestLBFGS_General(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x0 = self.system.get_random_minimized_configuration(tol=1).coords
        self.pot = self.system.get_potential()

    def test_event(self):
        self.called = False

        def event(coords=None, energy=None, rms=None):
            self.called = True

        opt = LBFGS(self.x0, self.pot, events=[event])
        opt.one_iteration()
        self.assertTrue(self.called)


class TestLBFGS_State(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()
        self.minimizer = LBFGS(self.x, self.pot)

    def test_state(self):
        # do several minimization iterations
        for i in range(10):
            self.minimizer.one_iteration()

        # get the state and save it
        ret = self.minimizer.get_result()
        state = self.minimizer.get_state()
        x1 = ret.coords.copy()

        # do several more iteration steps
        for i in range(10):
            self.minimizer.one_iteration()

        # now make a new minimizer and do several iterations
        minimizer2 = LBFGS(x1, self.pot)
        minimizer2.set_state(state)
        for i in range(10):
            minimizer2.one_iteration()

        # test that the two minimizers are in the same state
        ret1 = self.minimizer.get_result()
        ret2 = minimizer2.get_result()
        self.assertEqual(ret1.energy, ret2.energy)
        self.assertTrue((ret1.coords == ret2.coords).all())

        state1 = self.minimizer.get_state()
        state2 = minimizer2.get_state()

        self.assertTrue((state1.y == state2.y).all())
        self.assertTrue((state1.s == state2.s).all())
        self.assertTrue((state1.rho == state2.rho).all())
        self.assertTrue((state1.dXold == state2.dXold).all())
        self.assertTrue((state1.dGold == state2.dGold).all())
        self.assertEqual(state1.H0, state2.H0)
        self.assertEqual(state1.k, state2.k)

    def test_reset(self):
        # do several minimization iterations
        m1 = LBFGS(self.x, self.pot)
        for i in range(10):
            m1.one_iteration()

        # reset the minimizer and do it again
        m1.reset()
        e, g = self.pot.getEnergyGradient(self.x)
        m1.update_coords(self.x, e, g)
        for i in range(10):
            m1.one_iteration()

        # do the same number of steps of a new minimizer
        m2 = LBFGS(self.x, self.pot)
        for i in range(10):
            m2.one_iteration()

        # they should be the same (more or less)
        n = min(m1.k, m1.M)
        self.assertAlmostEqual(m1.H0, m2.H0, 5)
        self.assertEqual(m1.k, m2.k)
        arrays_nearly_equal(self, m1.y[:n, :], m2.y[:n, :])
        arrays_nearly_equal(self, m1.s[:n, :], m2.s[:n, :])
        arrays_nearly_equal(self, m1.rho[:n], m2.rho[:n])

        res1 = m1.get_result()
        res2 = m2.get_result()
        self.assertNotEqual(res1.nfev, res2.nfev)
        self.assertNotEqual(res1.nsteps, res2.nsteps)
        self.assertAlmostEqual(res1.energy, res2.energy)
        arrays_nearly_equal(self, res1.coords, res2.coords)


class TestLBFGS_discontinous(unittest.TestCase):
    def test1(self):
        pot = DiscontinuousHarmonic()
        x0 = np.array([-10, 1])
        opt = LBFGS(x0, pot, debug=True)
        res = opt.run()
        self.assertFalse(res.success)


class TestLBFGS_wolfe(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_minimized_configuration(tol=1e-1).coords
        self.pot = self.system.get_potential()

    def test(self):
        minimizer = LBFGS(self.x.copy(), self.pot, debug=True)
        minimizer._use_wolfe = True
        ret = minimizer.run()
        self.assertTrue(ret.success)

        print("\n\n")
        minimizer = LBFGS(self.x.copy(), self.pot, debug=True)
        ret_nowolfe = minimizer.run()
        self.assertTrue(ret_nowolfe.success)

        print("nfev wolfe, nowolfe", ret.nfev, ret_nowolfe.nfev, ret.energy, ret_nowolfe.energy)


class TestLBFGS_armijo(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()

    def test(self):
        minimizer = LBFGS(self.x.copy(), self.pot, armijo=True, debug=True)
        ret = minimizer.run()
        self.assertTrue(ret.success)

        print("\n\n")
        minimizer = LBFGS(self.x.copy(), self.pot, armijo=False, debug=True)
        ret_nowolfe = minimizer.run()
        self.assertTrue(ret_nowolfe.success)

        self.assertAlmostEqual(ret.energy, ret_nowolfe.energy, delta=1e-3)

        print("nfev armijo, noarmijo", ret.nfev, ret_nowolfe.nfev, ret.energy, ret_nowolfe.energy)


class TestLBFGSCython(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()

    def test(self):
        minimizer = LBFGS(self.x.copy(), self.pot, debug=True)
        minimizer._cython = True
        ret = minimizer.run()
        m2 = LBFGS(self.x.copy(), self.pot, debug=True)
        minimizer._cython = True
        ret2 = m2.run()

        print("cython", ret.nfev, ret2.nfev)
        self.assertEqual(ret.nfev, ret2.nfev)
        self.assertAlmostEqual(ret.energy, ret2.energy, 5)


class TestLBFGSFortran(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()

    def test(self):
        minimizer = LBFGS(self.x.copy(), self.pot, fortran=True, debug=True)
        ret = minimizer.run()
        m2 = LBFGS(self.x.copy(), self.pot, fortran=False, debug=True)
        ret2 = m2.run()

        print("fortran", ret.nfev, ret2.nfev)
        # self.assertEqual(ret.nfev, ret2.nfev)
        self.assertAlmostEqual(ret.energy, ret2.energy, 5)


if __name__ == "__main__":
    unittest.main()
        
        
        
