import unittest


import numpy as np

from pele.utils.vec3 import invert3x3
from pele.utils import rotations
from pele.utils.rotations import q2aa, mx2q, q_multiply, random_aa, random_q, q2mx

rot_epsilon = 1e-6


def aa2q(AA):
    """
    convert angle axis to quaternion
    input V: angle axis vector of lenth 3
    output Q: quaternion of length 4
    """
    q = np.zeros(4, np.float64)

    thetah = 0.5 * np.linalg.norm(AA)
    q[0] = np.cos(thetah)

    # do linear expansion for small epsilon
    if thetah < rot_epsilon:
        q[1:] = 0.5 * AA
    else:
        q[1:] = 0.5 * np.sin(thetah) * AA / thetah

    # make sure to have normal form
    if q[0] < 0.0: q = -q
    return q


def mx2aa(m):
    return q2aa(mx2q(m))


def aa2mx(p):
    return q2mx(aa2q(p))


def rotate_aa(p1, p2):
    """
    change a given angle axis rotation p2 by the
    rotation p1
    """
    return q2aa(q_multiply(aa2q(p2), aa2q(p1)))


class TestRotations(unittest.TestCase):
    def arrays_equal(self, a1, a2, places=5):
        self.assertEqual(a1.shape, a2.shape)
        for v1, v2 in zip(a1.reshape(-1), a2.reshape(-1)):
            self.assertAlmostEqual(v1, v2, places=places)

    def test_aa2q(self):
        aa = random_aa()
        q1 = aa2q(aa)
        q2 = rotations.aa2q(aa)
        self.arrays_equal(q1, q2)

    def test_rotate_aa(self):
        p1 = random_aa()
        p2 = random_aa()
        p3 = rotate_aa(p1, p2)
        p4 = rotations.rotate_aa(p1, p2)
        self.arrays_equal(p3, p4)

    def test_mx2aa(self):
        mx = rotations.q2mx(random_q())
        p1 = mx2aa(mx)
        p2 = rotations.mx2aa(mx)
        self.arrays_equal(p1, p2)

    def test_aa2mx(self):
        aa = random_aa()
        mx1 = aa2mx(aa)
        mx2 = rotations.aa2mx(aa)
        self.arrays_equal(mx1, mx2)

    def test_mx2q_fork1(self):
        mx = np.array([[0.84678162, 0.50684295, 0.16146551],
                       [-0.01866944, -0.27503628, 0.96125257],
                       [0.53161296, -0.81698548, -0.2234332]])
        qtrue = -np.array([-0.58058422, 0.76571063, 0.15938577, 0.22628603])
        q = rotations.mx2q(mx)
        self.arrays_equal(q, qtrue, 3)

    def test_mx2q_fork2(self):
        mx = np.array([[0.60801609, -0.66494404, 0.43378088],
                       [0.02750384, -0.52840591, -0.84854625],
                       [0.79344815, 0.52786041, -0.30299078]])
        qtrue = -np.array([-0.44063006, -0.78093098, 0.20406419, -0.39287372])
        q = rotations.mx2q(mx)
        self.arrays_equal(q, qtrue, 3)

    def test_mx2q_fork3(self):
        mx = np.array([[-0.44019684, -0.78441865, -0.43693721],
                       [-0.70029416, -0.00463489, 0.71383934],
                       [-0.56197405, 0.6202144, -0.54728353]])
        qtrue = np.array([0.04439801, -0.52719104, 0.70406773, 0.4736951])
        q = rotations.mx2q(mx)
        self.arrays_equal(q, qtrue, 3)

    def test_mx2q_fork4(self):
        mx = np.array([[-0.16802893, -0.42028749, -0.89169765],
                       [0.71301072, -0.67644975, 0.18447614],
                       [-0.68072167, -0.60479265, 0.41333261]])
        qtrue = np.array([0.37711203, -0.52323231, -0.13986293, 0.75130075])
        q = rotations.mx2q(mx)
        self.arrays_equal(q, qtrue)


class TestVec3(unittest.TestCase):
    def test_invert3x3(self):
        q = rotations.random_q()
        mx = rotations.q2mx(q)
        mxi1 = invert3x3(mx)
        mxi2 = np.linalg.inv(mx)
        self.assertEqual(mxi1.shape, mxi2.shape)
        for v1, v2 in zip(mxi1.reshape(-1), mxi2.reshape(-1)):
            self.assertAlmostEqual(v1, v2, places=5)


if __name__ == "__main__":
    unittest.main()

