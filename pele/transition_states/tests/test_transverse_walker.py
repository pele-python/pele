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
        self.assertLess(np.dot(dx, self.evec), 1e-6)


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

