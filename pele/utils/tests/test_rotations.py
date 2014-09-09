import unittest
from itertools import izip

import numpy as np

from pele.utils import rotations
from pele.utils.rotations import q2aa, mx2q, q_multiply, random_aa, random_q

rot_epsilon = 1e-6

def aa2q( AA ):
    """
    convert angle axis to quaternion
    input V: angle axis vector of lenth 3
    output Q: quaternion of length 4
    """
    q = np.zeros(4, np.float64)

    thetah = 0.5 * np.linalg.norm( AA ) 
    q[0]  = np.cos( thetah )

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

def rotate_aa(p1, p2):
    """
    change a given angle axis rotation p1 by the
    rotation p2
    """
    return q2aa(q_multiply( aa2q(p2), aa2q(p1) ))


class TestRotations(unittest.TestCase):
    def arrays_equal(self, a1, a2, places=5):
        self.assertEqual(a1.shape, a2.shape)
        for v1, v2 in izip(a1.reshape(-1), a2.reshape(-1)):
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
        

if __name__ == "__main__":
    unittest.main()
