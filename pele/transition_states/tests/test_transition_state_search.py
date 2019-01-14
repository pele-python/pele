from __future__ import print_function
import unittest


import numpy as np
import matplotlib.pyplot as plt

from pele.transition_states import FindTransitionState, findTransitionState
from pele.potentials import BasePotential
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
    x, y = np.meshgrid(np.arange(-2, 2, .1), np.arange(-2, 2, .1))
    pot = SimpleTSPot()
    energies = [pot.getEnergy(np.array([xi, yi])) for xi, yi in zip(x.reshape(-1), y.reshape(-1))]
    energies = np.array(energies).reshape(x.shape)
    plt.contourf(x, y, energies)
    plt.show()


class TestSimpleTSPot(_base_test._TestConfiguration):
    def setUp(self):
        self.pot = SimpleTSPot()
        self.x0 = np.array([.1, 1])
        self.e0 = -0.27335478531821539


class TestHarmonicPot(_base_test._TestConfiguration):
    def setUp(self):
        self.pot = HarmonicPot()
        self.x0 = np.array([.1, 1])
        self.e0 = 2.11


def print_event(coords=None, **kwargs):
    print("coords", coords)


class TestFindTransitionStateSimplePot(unittest.TestCase):
    def setUp(self):
        self.pot = SimpleTSPot()
        self.x0 = np.array([.1, .1])
        self.xts = np.zeros(self.x0.size)
        self.ets = self.pot.getEnergy(self.xts)

    def test1(self):
        # plot_pot()
        opt = FindTransitionState(self.x0, self.pot, orthogZeroEigs=None,
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
        ret = findTransitionState(self.x0, self.pot, orthogZeroEigs=None)
        self.assertTrue(ret.success)
        assert_arrays_almost_equal(self, ret.coords, self.xts, places=3)

    def test_2(self):
        self.called = False

        def event(**kwargs):
            self.called = True

        opt = FindTransitionState(self.x0, self.pot, orthogZeroEigs=None,
                                  tangentSpaceQuenchParams=dict(maxstep=1.),
                                  event=event)
        ret = opt.run()
        self.assertTrue(ret.success)
        self.assertTrue(self.called)

    def test_from_near_minimum(self):
        print("\n\nstarting from a minimum")
        x0 = np.array([.6, .1])
        opt = FindTransitionState(x0, self.pot, orthogZeroEigs=None,
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
        x0 = np.array([.6, .1])
        opt = FindTransitionState(x0, self.pot, orthogZeroEigs=None,
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
        x0 = np.array([.2, 0])
        opt = FindTransitionState(x0, pot, orthogZeroEigs=None,
                                  iprint=1,
                                  verbosity=10,  # event=print_event,
                                  hessian_diagonalization=True
        )
        ret = opt.run()
        self.assertFalse(ret.success)
        print(ret)


if __name__ == "__main__":
    unittest.main()

