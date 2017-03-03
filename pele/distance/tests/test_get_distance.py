from __future__ import division
import unittest
import numpy as np
from pele.distance import get_distance, Distance

def cround(r):
    if r > 0.0:
        r = np.floor(r + 0.5)
    else:
        r = np.ceil(r - 0.5)
    return r

class TestGetDistance(unittest.TestCase):
    """
    Test the get_distance method
    """

    # Method for manually calculating the distance in one dimension
    def _distance_1d(self, coord1, coord2, dim, use_leesedwards, shear):
        if use_leesedwards and dim in [0, 1]:
            d12 = [0, 0]
            d12[0] = coord1[0] - coord2[0]
            d12[1] = coord1[1] - coord2[1]

            round_y = cround(d12[1] / self.boxv[1])
            tmp12 = [d12[0] - round_y * shear * self.boxv[0],
                     d12[1] - round_y * self.boxv[1]]

            d12[0] -= cround(d12[0] / self.boxv[0]) * self.boxv[0]
            tmp12[0] -= cround(tmp12[0] / self.boxv[0]) * self.boxv[0]

            if d12[0] ** 2 + d12[1] ** 2 > tmp12[0] ** 2 + tmp12[1] ** 2:
                return tmp12[dim]
            else:
                return d12[dim]
        else:
            # Use distance to nearest image convention
            dist = coord1[dim] - coord2[dim]
            return dist - cround(dist / self.boxv[dim]) * self.boxv[dim]


    # Method for manually calculating the distance
    def _distance (self, coord1, coord2, ndim, use_leesedwards=False, shear=0.):
        return [self._distance_1d(coord1, coord2, dim, use_leesedwards, shear) for dim in xrange(ndim)]


    def setUp(self):
        boxlength = 10.
        self.test_repeat = 1000
        self.boxv = np.array([boxlength, boxlength, boxlength])
        self.x1s = np.random.uniform(-boxlength, boxlength, (self.test_repeat, 3))
        self.x2s = np.random.uniform(-boxlength, boxlength, (self.test_repeat, 3))


    def test_cartesian(self):
        dist_method = Distance.CARTESIAN
        for dim in [2, 3]:
            for x1, x2 in zip(self.x1s, self.x2s):
                diff = get_distance(x1, x2, dim, dist_method)
                for i in range(dim):
                    diff_target = x1[i] - x2[i]
                    self.assertAlmostEqual(diff[i], diff_target)


    def test_periodic(self):
        dist_method = Distance.PERIODIC
        for dim in [2, 3]:
            for x1, x2 in zip(self.x1s, self.x2s):
                diff = get_distance(x1, x2, dim, dist_method, self.boxv[:dim])
                diff_target = self._distance(x1, x2, dim)
                for i in range(dim):
                    self.assertAlmostEqual(diff[i], diff_target[i])


    def test_leesedwards(self):
        dist_method = Distance.LEES_EDWARDS
        shear = 0.1
        for dim in [2, 3]:
            for x1, x2 in zip(self.x1s, self.x2s):
                diff = get_distance(x1, x2, dim, dist_method, self.boxv[:dim], shear)
                diff_target = self._distance(x1, x2, dim, use_leesedwards=True, shear=shear)
                for i in range(dim):
                    self.assertAlmostEqual(diff[i], diff_target[i])


if __name__ == "__main__":
    unittest.main()
