from __future__ import print_function
import unittest


import numpy as np
#from numpy.testing import assert_array_almost_equal, assert_array_almost_equal_nulp

from pele.angleaxis._aadist import rmdrvt, sitedist_grad, sitedist
from pele.utils.rotations import vec_random
from pele.angleaxis._aa_utils import _rot_mat_derivative, _sitedist_grad, _sitedist
from pele.utils import rotations


class TestRmDrvt(unittest.TestCase):
    def assert_array_almost_equal(self, v1, v2, **kwargs):
        for x1, x2 in zip(v1.reshape(-1), v2.reshape(-1)):
            self.assertAlmostEqual(x1, x2, **kwargs)
        
    def test1(self):
        P = vec_random() * np.random.uniform(1e-5,1)
        self.check(P, True)

    def test_small_theta(self):
        P = vec_random() * np.random.uniform(0,1e-7)
        print(P)
        self.check(P, True)
    
    def check(self, P, with_grad):
        rm, drm1, drm2, drm3 = rmdrvt(P, with_grad)
        rmp, drm1p, drm2p, drm3p = _rot_mat_derivative(P, with_grad)
        print(rm, rmp)
        print("rm ", rm)
        print("rmp", rmp)
        print(np.abs(rm - rmp)/np.abs(rm))
        self.assert_array_almost_equal(rm, rmp)
        self.assert_array_almost_equal(drm1, drm1p, places=4)
        self.assert_array_almost_equal(drm2, drm2p, places=4)
        self.assert_array_almost_equal(drm3, drm3p, places=4)
    
    def test_big_small(self):
        P = vec_random()
        P1 = P * 1.1e-6
        P2 = P * 0.9e-6
#        print P1.dot(P1) > 1e-12, P2.dot(P2) > 1e-12

        rm, drm1, drm2, drm3 = rmdrvt(P1, True)
        rmp, drm1p, drm2p, drm3p = rmdrvt(P2, True)
#        print rm
#        print rmp
        self.assert_array_almost_equal(rm, rmp, places=4)
        self.assert_array_almost_equal(drm1, drm1p, places=4)
        self.assert_array_almost_equal(drm2, drm2p, places=4)
        self.assert_array_almost_equal(drm3, drm3p, places=4)


    def test_rot_mat1(self):
        p = np.array([1., 2, 3])
        p /= np.linalg.norm(p)
        print("test_rot_mat1")
        print(p)
        mx, r1, r2, r3 = _rot_mat_derivative(p, with_grad=True)
        print(mx)
        print(r1)
        print(r2)
        print(r3)

    def test_rot_mat_small_theta(self):
        p = np.array([1., 2, 3])
        p /= np.linalg.norm(p) * 1e7
        print("test_rot_mat1 small_theta")
        print(p)
        mx, r1, r2, r3 = _rot_mat_derivative(p, with_grad=True)
        print(mx)
        print(r1)
        print(r2)
        print(r3)

class TestSiteDistGrad(unittest.TestCase):
    def assert_array_almost_equal(self, v1, v2, **kwargs):
        for x1, x2 in zip(v1.reshape(-1), v2.reshape(-1)):
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

class TestRotations(unittest.TestCase):
    def test_mx2q(self):
        mx = np.array(list(range(9))).reshape([3,3])
        qnew = rotations.mx2q(mx)
        print(repr(qnew))
        qtrue = np.array([ 1.80277564,  0.2773501 , -0.5547002 ,  0.2773501 ])
        for v1, v2 in zip(qnew, qtrue):
            self.assertAlmostEqual(v1, v2, 4)
    
    def test_aa2q(self):
        print("\ntest_aa2q")
        p = np.array(list(range(1,4)), dtype=float)
        print(p)
        p /= np.linalg.norm(p)

        q = rotations.aa2q(p)
        print(repr(q))
        qtrue = np.array([ 0.87758256,  0.12813186,  0.25626373,  0.38439559])
        for v1, v2 in zip(q, qtrue):
            self.assertAlmostEqual(v1, v2, 4)

        
    def test_q2aa(self):
        print("\ntest_q2aa")
        v = np.array(list(range(1,4)), dtype=float)
        v /= np.linalg.norm(v)
        q = np.zeros(4)
        q[1:4] = v
        q[0] = 4
        print(q)
        aa = rotations.q2aa(q)
        print(repr(aa))
        aatrue = np.array([ 0.1309466 ,  0.26189321,  0.39283981])
        for v1, v2 in zip(aa, aatrue):
            self.assertAlmostEqual(v1, v2, 4)
        
    def test_mx2aa(self):
        print("test mx2aa")
        mx = np.array(list(range(9))).reshape([3,3])
        aa = rotations.mx2aa(mx)
        print(repr(aa))
        aatrue = np.array([ 0.29425463, -0.58850926,  0.29425463])
        for v1, v2 in zip(aa, aatrue):
            self.assertAlmostEqual(v1, v2, 4)
 
    def test_q_multiply(self):
        print("\ntest q_multiply")
        q1 = np.array(list(range(1,5)), dtype=float)
        q2 = np.array(list(range(2,6)), dtype=float)
        print(q1)
        q3 = rotations.q_multiply(q1, q2)
        print(repr(q3))
        qtrue = np.array([-36.,   6.,  12.,  12.])
        for v1, v2 in zip(q3, qtrue):
            self.assertAlmostEqual(v1, v2, 4)
    
    def test_rotate_aa(self):
        print("\ntest rotate_aa")
        p1 = np.array(list(range(1,4)), dtype=float)
        p2 = p1 + 1
        p3 = rotations.rotate_aa(p1, p2)
        print(repr(p3))
        ptrue = np.array([ 0.74050324,  1.64950785,  2.20282887])
        for v1, v2 in zip(p3, ptrue):
            self.assertAlmostEqual(v1, v2, 4)

        
if __name__ == "__main__":
    unittest.main()

