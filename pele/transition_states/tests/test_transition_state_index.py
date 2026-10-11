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
            orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertFalse(ret.success)

    def test_accepts_index_one_when_zero_modes_leave_no_transverse_direction(self):
        def remove_zero_mode(vector, coords):
            return np.array([vector[0], 0.0])

        ret = FindTransitionState(
            np.zeros(2), _QuadraticPotential([[-2.0, 0.0], [0.0, 0.0]]),
            eigenvec0=np.array([1.0, 0.0]), orthogZeroEigs=remove_zero_mode,
            nsteps=1, verbosity=0,
        ).run()
        self.assertTrue(ret.success)
        self.assertLess(ret.rms, 1e-4)

    def test_accepts_one_dimensional_index_one_stationary_point(self):
        ret = FindTransitionState(
            np.zeros(1), _QuadraticPotential([[-2.0]]),
            eigenvec0=np.ones(1), orthogZeroEigs=None, nsteps=1, verbosity=0,
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
            orthogZeroEigs=None, nsteps=1, verbosity=0,
        ).run()
        self.assertEqual(ret.nfev, pot.nfev)
        self.assertTrue(ret.success, msg=repr(ret))

    def test_success_requires_a_converged_returned_primary_eigenpair(self):
        search = FindTransitionState(
            np.zeros(1), _QuadraticPotential([[-2.0]]),
            eigenvec0=np.ones(1), orthogZeroEigs=None, nsteps=1, verbosity=0,
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
            orthogZeroEigs=None, demand_initial_negative_vec=True,
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
