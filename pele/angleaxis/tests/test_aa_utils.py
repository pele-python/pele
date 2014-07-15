import unittest
from itertools import izip

import numpy as np
#from numpy.testing import assert_array_almost_equal, assert_array_almost_equal_nulp

from pele.angleaxis._aadist import rmdrvt, sitedist_grad, sitedist
from pele.utils.rotations import vec_random
from pele.angleaxis._aa_utils import _rot_mat_derivative, _sitedist_grad, _sitedist


class TestRmDrvt(unittest.TestCase):
    def assert_array_almost_equal(self, v1, v2, **kwargs):
        for x1, x2 in izip(v1.reshape(-1), v2.reshape(-1)):
            self.assertAlmostEqual(x1, x2, **kwargs)
        
    def test1(self):
        P = vec_random() * np.random.uniform(1e-5,1)
        self.check(P, True)

    def test_small_theta(self):
        P = vec_random() * np.random.uniform(0,1e-7)
        print P
        self.check(P, True)
    
    def check(self, P, with_grad):
        rm, drm1, drm2, drm3 = rmdrvt(P, with_grad)
        rmp, drm1p, drm2p, drm3p = _rot_mat_derivative(P, with_grad)
        print rm, rmp
        print "rm ", rm
        print "rmp", rmp
        print np.abs(rm - rmp)/np.abs(rm)
        self.assert_array_almost_equal(rm, rmp)
        self.assert_array_almost_equal(drm1, drm1p, places=4)
        self.assert_array_almost_equal(drm2, drm2p, places=4)
        self.assert_array_almost_equal(drm3, drm3p, places=4)
    
    def test_big_small(self):
        P = vec_random()
        P1 = P * (1.1e-6)
        P2 = P * (0.9e-6)
#        print P1.dot(P1) > 1e-12, P2.dot(P2) > 1e-12

        rm, drm1, drm2, drm3 = rmdrvt(P1, True)
        rmp, drm1p, drm2p, drm3p = rmdrvt(P1, True)
#        print rm
#        print rmp
        self.assert_array_almost_equal(rm, rmp, places=4)
        self.assert_array_almost_equal(drm1, drm1p, places=4)
        self.assert_array_almost_equal(drm2, drm2p, places=4)
        self.assert_array_almost_equal(drm3, drm3p, places=4)


    def test_rot_mat1(self):
        p = np.array([1., 2, 3]);
        p /= np.linalg.norm(p);
        print "test_rot_mat1"
        print p;
        mx, r1, r2, r3 = _rot_mat_derivative(p, with_grad=True)
        print mx
        print r1
        print r2
        print r3

    def test_rot_mat_small_theta(self):
        p = np.array([1., 2, 3]);
        p /= np.linalg.norm(p) * 1e7;
        print "test_rot_mat1 small_theta"
        print p;
        mx, r1, r2, r3 = _rot_mat_derivative(p, with_grad=True)
        print mx
        print r1
        print r2
        print r3

class TestSiteDistGrad(unittest.TestCase):
    def assert_array_almost_equal(self, v1, v2, **kwargs):
        for x1, x2 in izip(v1.reshape(-1), v2.reshape(-1)):
            self.assertAlmostEqual(x1, x2, **kwargs)
        
    def test1(self):
        drij = np.random.uniform(-1,1,3)
        p1 = vec_random() * 0.5
        p2 = vec_random() * 0.5
        W = 1.3
        S = np.random.uniform(-1,1,[3,3])
        cog = np.random.uniform(-1,1,3)
        
        g_M, g_P = sitedist_grad(drij, p1, p2, S, W, cog)
        g_Mp, g_Pp = _sitedist_grad(drij, p1, p2, S, W, cog)
        self.assert_array_almost_equal(g_M, g_Mp, places=4)
        self.assert_array_almost_equal(g_P, g_Pp, places=4)

    def test2(self):
        # rotate around the z axis by pi/2
        P = np.array([0., 0., 1.]) * np.pi/2
        v1 = np.array([1.,0,0])
        rm = _rot_mat_derivative(P, False)[0]
#        rm = rmdrvt(P, False)[0]
        v2 = rm.dot(v1)
        v2_true = np.array([0.,1,0])
        self.assert_array_almost_equal(v2, v2_true)
    
    def test_sitedist(self):
        drij = np.random.uniform(-1,1,3)
        p1 = vec_random() * 0.5
        p2 = vec_random() * 0.5
        W = 1.3
        S = np.random.uniform(-1,1,[3,3])
        cog = np.random.uniform(-1,1,3)
        
        dist = sitedist(drij, p1, p2, S, W, cog)
        dist2 = _sitedist(drij, p1, p2, S, W, cog)
        self.assertAlmostEqual(dist, dist2, places=4)


if __name__ == "__main__":
    unittest.main()
