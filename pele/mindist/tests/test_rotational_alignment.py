import unittest

import numpy as np

from pele.mindist import (findrotation_kabsch, findrotation_kearsley, findrotation,
                          MeasureAtomicCluster, TransformAtomicCluster)
from pele.utils.rotations import random_q, q2mx, mx2q
from pele.potentials.tests._base_test import assert_arrays_almost_equal

def kabsch_wrap(x1, x2, align_com=True):
    """compute and return the distance to agree with findrotation_kearsley
    """
    mx = findrotation_kabsch(x1, x2, align_com=align_com)
    tform = TransformAtomicCluster()
    measure = MeasureAtomicCluster()
    x2 = x2.copy()
    x1 = x1.copy()
    tform.translate(x1, -measure.get_com(x1))
    tform.translate(x2, -measure.get_com(x2))
    tform.rotate(x2, mx)
    dist = measure.get_dist(x1, x2)
    return dist, mx
    
    

class TestFindRotation(unittest.TestCase):
    def setUp(self):
        natoms = 10
        self.x1 = np.random.uniform(0, 1, 3*natoms)
        self.x2 = np.random.uniform(0, 1, 3*natoms)
        
        self.measure = MeasureAtomicCluster(permlist=[])
        self.transform = TransformAtomicCluster(can_invert=True)
        
        self.findrot = findrotation
        
    def test1(self):
        rot = q2mx(random_q())
        x2 = self.x1.copy()
        self.transform.rotate(x2, rot)
        
        x1 = self.x1
        dist, mx = self.findrot(x1, x2)
        self.assertLess(dist, 1e-4)
        
        self.transform.translate(x2, -self.measure.get_com(x2))
        self.transform.rotate(x2, mx)
        self.transform.translate(x2, self.measure.get_com(x1))

        assert_arrays_almost_equal(self, x1, x2)


    def test2(self):
        rot = q2mx(random_q())
        dx = .01
        x2 = self.x1.copy()
        x2[0] += .01
        self.transform.rotate(x2, rot)
        
        x1old = self.x1.copy()
        x2old = x2.copy()
        dist, mx = self.findrot(self.x1, x2)
        # check it didn't change the arrays
        assert_arrays_almost_equal(self, x1old, self.x1)
        assert_arrays_almost_equal(self, x2old, x2)
        
        self.assertLess(dist, dx)
    
    def test_compare(self):
        dkear, rot1 = findrotation_kearsley(self.x1, self.x2, align_com=True)
        dkab, rot2 = kabsch_wrap(self.x1, self.x2, align_com=True)
        
        qkear = mx2q(rot1)
        qkab = mx2q(rot2)
        
        self.assertAlmostEqual(dkear, dkab, places=4)
        assert_arrays_almost_equal(self, qkear, qkab)
        
class TestFindRotKearsley(TestFindRotation):
    def setUp(self):
        TestFindRotation.setUp(self)
        self.findrot = findrotation_kearsley

class TestFindRotKabsch(TestFindRotation):
    def setUp(self):
        TestFindRotation.setUp(self)
        self.findrot = kabsch_wrap



if __name__ == "__main__":
    unittest.main()

