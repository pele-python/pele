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
        g = 2.0 * x
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
        rng_state = np.random.get_state()
        self.addCleanup(np.random.set_state, rng_state)
        np.random.seed(0)
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
        array_fields = ("s", "y", "rho", "dXold", "dGold")
        for field in array_fields:
            saved = getattr(state, field)
            original = getattr(self.minimizer, field)
            np.testing.assert_array_equal(saved, original)
            self.assertFalse(np.shares_memory(saved, original))
        self.assertEqual(
            (state.k, state.H0, state.have_dXold),
            (self.minimizer.k, self.minimizer.H0, self.minimizer._have_dXold),
        )

        # do several more iteration steps
        for i in range(10):
            self.minimizer.one_iteration()

        # now make a new minimizer and do several iterations
        minimizer2 = LBFGS(x1, self.pot)
        minimizer2.set_state(state)
        restored = minimizer2.get_state()
        for field in array_fields:
            np.testing.assert_array_equal(getattr(restored, field), getattr(state, field))
        self.assertEqual(
            (restored.k, restored.H0, restored.have_dXold),
            (state.k, state.H0, state.have_dXold),
        )
        for i in range(10):
            minimizer2.one_iteration()

        # Restored memory is exact; subsequent reductions can round differently.
        tolerance = dict(rtol=1e-12, atol=1e-12, equal_nan=False)
        ret1 = self.minimizer.get_result()
        ret2 = minimizer2.get_result()
        np.testing.assert_allclose(ret1.energy, ret2.energy, **tolerance)
        np.testing.assert_allclose(ret1.coords, ret2.coords, **tolerance)
        np.testing.assert_allclose(ret1.grad, ret2.grad, **tolerance)

        state1 = self.minimizer.get_state()
        state2 = minimizer2.get_state()
        for field in array_fields:
            np.testing.assert_allclose(getattr(state1, field), getattr(state2, field), **tolerance)
        np.testing.assert_allclose(state1.H0, state2.H0, **tolerance)
        self.assertEqual(state1.k, state2.k)
        self.assertEqual(state1.have_dXold, state2.have_dXold)

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

        print(
            "nfev wolfe, nowolfe",
            ret.nfev,
            ret_nowolfe.nfev,
            ret.energy,
            ret_nowolfe.energy,
        )


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

        print(
            "nfev armijo, noarmijo",
            ret.nfev,
            ret_nowolfe.nfev,
            ret.energy,
            ret_nowolfe.energy,
        )


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
