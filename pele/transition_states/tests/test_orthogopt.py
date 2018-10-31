import unittest

import numpy as np

from pele.transition_states._orthogopt import orthogopt, orthogopt_slow, orthogopt_translation_only
from pele.transition_states._zeroev import zeroEV_translation, zeroEV_rotation


class TestOrthogopt(unittest.TestCase):
    def test_o(self):
        # test orthogopt
        from pele.transition_states import zeroEV_cluster

        natoms = 13
        vec = np.random.uniform(-1, 1, natoms * 3)
        vec /= np.linalg.norm(vec)
        coords = np.random.uniform(-1, 1, natoms * 3)

        # orthogonalize vec
        vec = orthogopt(vec, coords)

        zeroev = zeroEV_cluster(coords)
        for u in zeroev:
            # print np.dot(u, vec)
            self.assertAlmostEqual(0., np.dot(u, vec), 5)

    def test_slow(self):
        # test orthogopt_slow
        from pele.transition_states import zeroEV_cluster

        natoms = 13
        vec = np.random.uniform(-1, 1, natoms * 3)
        vec /= np.linalg.norm(vec)
        coords = np.random.uniform(-1, 1, natoms * 3)

        # orthogonalize vec
        vec = orthogopt_slow(vec, coords)

        zeroev = zeroEV_cluster(coords)
        for u in zeroev:
            # print np.dot(u, vec)
            self.assertAlmostEqual(0., np.dot(u, vec), 5)

    def test_translation_only(self):
        # test orthogopt with translations only
        natoms = 13
        vec = np.random.uniform(-1, 1, natoms * 3)
        vec /= np.linalg.norm(vec)
        coords = np.random.uniform(-1, 1, natoms * 3)

        # orthogonalize vec
        vec = orthogopt_translation_only(vec, coords)

        zeroev = zeroEV_translation(coords)
        for u in zeroev:
            # print np.dot(u, vec)
            self.assertAlmostEqual(0., np.dot(u, vec), 5)

        rotev = zeroEV_rotation(coords)
        for u in rotev:
            # print np.dot(u, vec)
            self.assertNotAlmostEqual(0., np.dot(u, vec), 5)


if __name__ == "__main__":
    unittest.main()
