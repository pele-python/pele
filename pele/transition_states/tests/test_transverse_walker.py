import unittest
import numpy as np

from pele.systems import LJCluster
from pele.transition_states._transverse_walker import _TransverseWalker
from pele.utils.rotations import vec_random_ndim


class TestTransverseWalker(unittest.TestCase):
    def setUp(self):
        natoms = 18
        self.system = LJCluster(natoms)
        self.pot = self.system.get_potential()
        self.evec = vec_random_ndim(3 * natoms)
        self.evec /= np.linalg.norm(self.evec)
        self.x = self.system.get_random_configuration()

    def test(self):
        """assert the motion is perpendicular to evec"""
        x0 = self.x.copy()
        opt = _TransverseWalker(self.x, self.pot, self.evec)
        ret = opt.run(10)
        xnew = ret.coords

        dx = xnew - x0
        self.assertLess(abs(np.dot(dx, self.evec)), 1e-6)


class TestTransverseWalkerHistory(unittest.TestCase):
    def test_unchanged_plane_keeps_the_same_trajectory(self):
        class QuadraticPotential:
            def getEnergyGradient(self, coords):
                gradient = np.array([1.0, 10.0, 100.0]) * coords
                return 0.5 * np.dot(coords, gradient), gradient

        pot = QuadraticPotential()
        coords = np.array([1.0, -0.7, 0.4])
        axis = np.array([0.0, 0.0, 1.0])
        baseline = _TransverseWalker(coords, pot, axis, tol=1e-12)
        expected = baseline.run(4).coords.copy()
        for sign in (1.0, -1.0):
            with self.subTest(sign=sign):
                opt = _TransverseWalker(coords, pot, axis, tol=1e-12)
                current = opt.run(3).coords.copy()
                energy, gradient = pot.getEnergyGradient(current)
                opt.update_eigenvec(sign * axis, -1.0)
                opt.update_coords(current, energy, gradient)
                np.testing.assert_allclose(opt.run(1).coords, expected,
                                           rtol=0.0, atol=1e-14)

    def test_rotation_restores_configured_initial_scale(self):
        class QuadraticPotential:
            def getEnergyGradient(self, coords):
                gradient = np.array([1.0, 2.0, 3.0]) * coords
                return 0.5 * np.dot(coords, gradient), gradient

        pot = QuadraticPotential()
        opt = _TransverseWalker(
            np.array([1.0, -0.7, 0.4]), pot, np.array([0.0, 0.0, 1.0]),
            tol=1e-12, H0=0.7, maxstep=10.0,
        )
        opt.run(3)
        axis = np.array([1.0, 0.0, 0.0])
        coords = np.array([-0.5, 1.0, 1.0])
        energy, gradient = pot.getEnergyGradient(coords)
        opt.update_eigenvec(axis, -1.0)
        opt.update_coords(coords, energy, gradient)
        fresh = _TransverseWalker(coords, pot, axis, tol=1e-12,
                                  H0=0.7, maxstep=10.0)
        np.testing.assert_allclose(opt.run(1).coords, fresh.run(1).coords,
                                   rtol=0.0, atol=1e-14)

    def test_motion_stays_in_updated_plane(self):
        class QuadraticPotential:
            hessian = np.array([[3.0, 1.0, 0.5],
                                [1.0, 2.0, 0.25],
                                [0.5, 0.25, 1.0]])

            def getEnergyGradient(self, coords):
                gradient = self.hessian @ coords
                return 0.5 * np.dot(coords, gradient), gradient

        pot = QuadraticPotential()
        opt = _TransverseWalker(
            np.array([1.0, -0.7, 0.4]), pot, np.array([0.0, 0.0, 1.0]),
            tol=1e-12,
        )
        # Build valid secant history in the original transverse plane.
        opt.run(3)
        axis = np.array([1.0, 0.0, 0.0])
        coords = np.array([-0.5, 1.0, 1.0])
        energy, gradient = pot.getEnergyGradient(coords)
        opt.update_eigenvec(axis, -1.0)
        opt.update_coords(coords, energy, gradient)

        start = coords.copy()
        for _ in range(5):
            previous = coords.copy()
            coords = opt.run(1).coords.copy()
            self.assertLess(abs(np.dot(coords - previous, axis)), 1e-14)
        self.assertGreater(np.linalg.norm(coords - start), 1e-5)


class TestTransverseWalker_NFEV(unittest.TestCase):
    def setUp(self):
        from pele.optimize.tests.test_nfev import _PotWrapper

        self.system = LJCluster(18)
        self.pot = _PotWrapper(self.system.get_potential())

    def test(self):
        x = self.system.get_random_configuration()
        evec = vec_random_ndim(x.size)
        evec /= np.linalg.norm(evec)
        opt = _TransverseWalker(x, self.pot, evec)
        ret = opt.run(10)
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertGreater(ret.nfev, 0)


if __name__ == "__main__":
    unittest.main()
