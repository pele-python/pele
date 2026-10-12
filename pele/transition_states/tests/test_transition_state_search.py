import unittest
from pathlib import Path


import numpy as np
import matplotlib.pyplot as plt

from pele.transition_states import FindTransitionState, findTransitionState
from pele.potentials import BasePotential, LJ
from pele.utils.xyz import read_xyz
from pele.potentials.tests import _base_test
from pele.potentials.tests._base_test import assert_arrays_almost_equal


class SimpleTSPot(BasePotential):
    nfev = 0
    # def __init__(self):
    # self.nfcalls = 0
    def getEnergyGradient(self, x):
        self.nfev += 1
        grad = np.zeros(x.size)

        dx = x.copy()
        dx[0] += 1
        e1 = -np.exp(-np.dot(dx, dx))
        grad -= 2 * dx * e1

        dx = x.copy()
        dx[0] -= 1
        e2 = -np.exp(-np.dot(dx, dx))
        grad -= 2 * dx * e2

        return e1 + e2, grad

    def getEnergy(self, x):
        return self.getEnergyGradient(x)[0]


class HarmonicPot(BasePotential):
    def getEnergy(self, x):
        return np.dot(x, x) + x.sum()

    def getEnergyGradient(self, x):
        e = self.getEnergy(x)
        return e, 2 * x + 1


def plot_pot():
    x, y = np.meshgrid(np.arange(-2, 2, 0.1), np.arange(-2, 2, 0.1))
    pot = SimpleTSPot()
    energies = [
        pot.getEnergy(np.array([xi, yi]))
        for xi, yi in zip(x.reshape(-1), y.reshape(-1))
    ]
    energies = np.array(energies).reshape(x.shape)
    plt.contourf(x, y, energies)
    plt.show()


class TestSimpleTSPot(_base_test.BaseTestCases._TestConfiguration):
    def setUp(self):
        self.pot = SimpleTSPot()
        self.x0 = np.array([0.1, 1])
        self.e0 = -0.27335478531821539


class TestHarmonicPot(_base_test.BaseTestCases._TestConfiguration):
    def setUp(self):
        self.pot = HarmonicPot()
        self.x0 = np.array([0.1, 1])
        self.e0 = 2.11


def print_event(coords=None, **kwargs):
    print("coords", coords)


class TestFindTransitionStateSimplePot(unittest.TestCase):
    def setUp(self):
        self.addCleanup(np.random.set_state, np.random.get_state())
        np.random.seed(0)
        self.pot = SimpleTSPot()
        self.x0 = np.array([0.1, 0.1])
        self.xts = np.zeros(self.x0.size)
        self.ets = self.pot.getEnergy(self.xts)

    def test1(self):
        # plot_pot()
        opt = FindTransitionState(
            self.x0,
            self.pot,
            index_seed=0, orthogZeroEigs=None,
            # iprint=1,
            # verbosity=10, event=print_event,
            # tol=1e-3,
            # lowestEigenvectorQuenchParams=dict(iprint=1, events=[print_event])
        )
        ret = opt.run()
        self.assertTrue(ret.success)
        assert_arrays_almost_equal(self, ret.coords, self.xts, places=3)
        self.assertAlmostEqual(ret.energy, self.ets, delta=1e-3)
        self.assertLess(ret.rms, 1e-3)
        self.assertEqual(ret.nfev + 1, self.pot.nfev)

    def test_wrapper(self):
        ret = findTransitionState(self.x0, self.pot, index_seed=0, orthogZeroEigs=None)
        self.assertTrue(ret.success)
        assert_arrays_almost_equal(self, ret.coords, self.xts, places=3)

    def test_2(self):
        self.called = False

        def event(**kwargs):
            self.called = True

        opt = FindTransitionState(
            self.x0,
            self.pot,
            index_seed=0, orthogZeroEigs=None,
            tangentSpaceQuenchParams=dict(maxstep=1.0),
            event=event,
        )
        ret = opt.run()
        self.assertTrue(ret.success)
        self.assertTrue(self.called)

    def test_from_near_minimum(self):
        print("\n\nstarting from a minimum")
        x0 = np.array([0.6, 0.1])
        opt = FindTransitionState(
            x0,
            self.pot,
            index_seed=0, orthogZeroEigs=None,
            iprint=1,
            verbosity=10,  # event=print_event,
            # tol=1e-3,
            # lowestEigenvectorQuenchParams=dict(iprint=1, events=[print_event])
        )
        ret = opt.run()
        print(ret)
        self.assertTrue(ret.success)
        assert_arrays_almost_equal(self, ret.coords, self.xts, places=3)

    def test_from_near_minimum_demand_negative_eigenvalue(self):
        print("\n\nstarting from a minimum demand")
        # demand that the eigenvalue is negative initially.
        # this should fail right away
        x0 = np.array([0.6, 0.1])
        opt = FindTransitionState(
            x0,
            self.pot,
            index_seed=0, orthogZeroEigs=None,
            demand_initial_negative_vec=True,
            iprint=1,
            verbosity=10,  # event=print_event,
            # tol=1e-3,
            # lowestEigenvectorQuenchParams=dict(iprint=1, events=[print_event])
        )
        ret = opt.run()
        print(ret)
        self.assertFalse(ret.success)
        self.assertEqual(ret.nsteps, 0)


class TestFindTS_BadPotential(unittest.TestCase):
    def test1(self):
        print("\n\ntesting find ts with harmonic potential")
        pot = HarmonicPot()
        x0 = np.array([0.2, 0])
        opt = FindTransitionState(
            x0,
            pot,
            index_seed=0, orthogZeroEigs=None,
            iprint=1,
            verbosity=10,  # event=print_event,
            hessian_diagonalization=True,
        )
        ret = opt.run()
        self.assertFalse(ret.success)
        print(ret)


class IndexTwoPotential(BasePotential):
    def getEnergyGradient(self, coords):
        curvature = np.array([-2.0, -1.0, 1.0])
        return 0.5 * np.dot(curvature * coords, coords), curvature * coords

    def getHessian(self, coords):
        raise AssertionError("Default index validation must not form a dense Hessian")


class RecoverableIndexTwoPotential(BasePotential):
    def getEnergyGradient(self, coords):
        x, y, z = coords
        return (-x*x + 0.25*(y*y - 1.0)**2 + 0.5*z*z,
                np.array([-2.0*x, y*(y*y - 1.0), z]))


class TestTransitionStateIndex(unittest.TestCase):
    def test_rejects_index_two_stationary_point(self):
        ret = FindTransitionState(
            np.zeros(3), IndexTwoPotential(), eigenvec0=np.array([1., 0., 0.]),
            index_seed=0, orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertFalse(ret.success)

    def test_rejects_recorded_lj13_index_two_saddle(self):
        with Path(__file__).with_name("lj13_index_two.xyz").open() as source:
            coords = read_xyz(source).coords.ravel()
        pot = LJ()
        values, vectors = np.linalg.eigh(pot.getHessian(coords))
        self.assertLess(values[1], -1.0)
        ret = FindTransitionState(
            coords, pot, eigenvec0=vectors[:, 0], index_seed=0, nsteps=1, verbosity=0,
        ).run()
        self.assertFalse(ret.success)

    def test_escapes_extra_unstable_mode(self):
        ret = FindTransitionState(
            np.zeros(3), RecoverableIndexTwoPotential(),
            eigenvec0=np.array([1., 0., 0.]), index_seed=0, orthogZeroEigs=None, verbosity=0,
        ).run()
        self.assertTrue(ret.success)
        self.assertAlmostEqual(abs(ret.coords[1]), 1.0, delta=1e-3)
        self.assertLess(ret.rms, 1e-4)


if __name__ == "__main__":
    unittest.main()
