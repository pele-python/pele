import unittest
from unittest import mock

import numpy as np

from pele.potentials import BasePotential
from pele.transition_states import FindTransitionState


class _QuadraticPotential(BasePotential):
    def __init__(self, hessian):
        self.hessian = np.asarray(hessian, dtype=float)
        self.nfev = 0

    def getEnergyGradient(self, coords):
        self.nfev += 1
        gradient = self.hessian @ coords
        return 0.5 * np.dot(coords, gradient), gradient

    def getHessian(self, coords):
        raise AssertionError("Default index validation must not form a dense Hessian")


class TestTransitionStateIndexEdges(unittest.TestCase):
    def test_rejects_rotated_index_two_with_stable_initial_rr_direction(self):
        # Eigenvalues are -2, -1, +1; the first transverse probe points along +1.
        pot = _QuadraticPotential([
            [-2.0, 0.0, 0.0],
            [0.0, -0.918372375740774, -0.3957172974185572],
            [0.0, -0.3957172974185572, 0.918372375740774],
        ])
        ret = FindTransitionState(
            np.zeros(3), pot, eigenvec0=np.array([1.0, 0.0, 0.0]),
            index_seed=0, orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertFalse(ret.success)

    def test_random_start_finds_index_two_mode_missed_by_fixed_guesses(self):
        # Spectrum -2, -1, +1, +1; both old fixed probes miss the -1 mode.
        pot = _QuadraticPotential([
            [-2.0, 0.0, 0.0, 0.0],
            [0.0, -0.07552173356870506, -0.3733085508688484, 0.9246281380138601],
            [0.0, -0.3733085508688484, 0.870426352344002, 0.3209340913541529],
            [0.0, 0.9246281380138601, 0.3209340913541529, 0.20509538122470305],
        ])
        np.testing.assert_allclose(np.linalg.eigvalsh(pot.hessian),
                                   [-2.0, -1.0, 1.0, 1.0], atol=1e-14)
        ret = FindTransitionState(
            np.zeros(4), pot, eigenvec0=np.array([1.0, 0.0, 0.0, 0.0]),
            index_seed=1, orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertEqual(ret.nfev, pot.nfev)
        self.assertFalse(ret.success)

    def test_weak_high_dimensional_mode_is_not_lost_to_rms_stopping(self):
        dimension = 4096
        hidden = np.random.default_rng(314159).normal(size=dimension)
        hidden[0] = 0.0
        hidden /= np.linalg.norm(hidden)

        class RankTwoPotential(BasePotential):
            def getEnergyGradient(self, coords):
                # Spectrum -2, -0.01, +0.01... without allocating a Hessian.
                gradient = 0.01 * coords - 0.02 * np.dot(hidden, coords) * hidden
                gradient[0] -= 2.01 * coords[0]
                return 0.5 * np.dot(coords, gradient), gradient

        primary = np.zeros(dimension)
        primary[0] = 1.0
        search = FindTransitionState(
            np.zeros(dimension), RankTwoPotential(), eigenvec0=primary,
            index_seed=0, orthogZeroEigs=None, verbosity=0,
        )
        result = search._get_lowest_transverse_eigenvector(
            np.zeros(dimension), np.zeros(dimension),
        )
        self.assertTrue(result.success)
        self.assertLess(result.eigenval, 0.0)

    def test_accepts_index_one_when_zero_modes_leave_no_transverse_direction(self):
        def remove_zero_mode(vector, coords):
            return np.array([vector[0], 0.0])

        ret = FindTransitionState(
            np.zeros(2), _QuadraticPotential([[-2.0, 0.0], [0.0, 0.0]]),
            eigenvec0=np.array([1.0, 0.0]), index_seed=0, orthogZeroEigs=remove_zero_mode,
            nsteps=1, verbosity=0,
        ).run()
        self.assertTrue(ret.success)
        self.assertLess(ret.rms, 1e-4)

    def test_accepts_one_dimensional_index_one_stationary_point(self):
        ret = FindTransitionState(
            np.zeros(1), _QuadraticPotential([[-2.0]]),
            eigenvec0=np.ones(1), index_seed=0, orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertTrue(ret.success)
        self.assertAlmostEqual(ret.eigenval, -2.0)
        self.assertLess(ret.rms, 1e-4)

    def test_accepts_stiff_index_one_stationary_point(self):
        # One negative mode and eleven positive modes with disparate stiffnesses.
        curvatures = np.r_[-2.0, np.geomspace(1.0, 1000.0, 11)]
        pot = _QuadraticPotential(np.diag(curvatures))
        primary = np.zeros(curvatures.size)
        primary[0] = 1.0
        ret = FindTransitionState(
            np.zeros(curvatures.size), pot, eigenvec0=primary,
            index_seed=0, orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertEqual(ret.nfev, pot.nfev)
        self.assertTrue(ret.success, msg=repr(ret))

    def test_success_requires_a_converged_returned_primary_eigenpair(self):
        search = FindTransitionState(
            np.zeros(1), _QuadraticPotential([[-2.0]]),
            eigenvec0=np.ones(1), index_seed=0, orthogZeroEigs=None, nsteps=1, verbosity=0,
        )
        find_primary = search._get_lowest_eigenvector_RR
        calls = 0

        def fail_redundant_refresh(*args, **kwargs):
            nonlocal calls
            result = find_primary(*args, **kwargs)
            calls += 1
            if calls >= 3:
                result.success = False
            return result

        with mock.patch.object(search, "_get_lowest_eigenvector_RR",
                               side_effect=fail_redundant_refresh):
            ret = search.run()
        # Skipping the redundant refresh is valid; accepting its failed result is not.
        self.assertTrue(not ret.success or search.leig_result.success)

    def test_recovery_failure_returns_energy_and_gradient_at_returned_coords(self):
        pot = _QuadraticPotential(np.diag([-2.0, -1.0, 1.0]))
        search = FindTransitionState(
            np.zeros(3), pot, eigenvec0=np.array([1.0, 0.0, 0.0]),
            index_seed=0, orthogZeroEigs=None, demand_initial_negative_vec=True,
            nsteps=2, verbosity=0,
        )
        search.npositive_max = 0
        find_primary = search._getLowestEigenVector

        def abort_after_recovery(coords, iteration, gradient=None):
            overlap = find_primary(coords, iteration, gradient=gradient)
            if iteration > 0:
                search.eigenval = 1.0
            return overlap

        with mock.patch.object(search, "_getLowestEigenVector",
                               side_effect=abort_after_recovery):
            ret = search.run()
        self.assertFalse(ret.success)
        self.assertGreater(np.linalg.norm(ret.coords), 1e-5)
        self.assertEqual(ret.nfev, pot.nfev)
        energy, gradient = pot.getEnergyGradient(ret.coords)
        self.assertAlmostEqual(ret.energy, energy, places=12)
        np.testing.assert_allclose(ret.grad, gradient, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(ret.rms, np.linalg.norm(gradient) / np.sqrt(3.0), places=12)


if __name__ == "__main__":
    unittest.main()
