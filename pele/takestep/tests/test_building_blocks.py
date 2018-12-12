import unittest

import numpy as np

from pele.utils import rotations
from pele.takestep.buildingblocks import (reduced_coordinates_displace, rotate, 
                                          uniform_displace)

def relative_angle(aa1, aa2):
    m1 = rotations.aa2mx(aa1)
    m2 = rotations.aa2mx(aa2)
    m2 = m1.transpose().dot(m2)
    theta = np.linalg.norm(rotations.mx2aa(m2))
    return theta
    

class TestUniformDisplace(unittest.TestCase):
    def do(self, indices=None):
        natoms = 4
        x = np.zeros(natoms*3)
        
        stepsize=1.5
        uniform_displace(stepsize, x, indices=indices)
        
        x = x.reshape(-1,3)
        rlist = [np.linalg.norm(r) for r in x]
        for r in rlist:
            self.assertLessEqual(r, stepsize)
        
        if indices is not None:
            for i, r in enumerate(rlist):
                if i in indices:
                    self.assertGreater(r, 0)
                else:
                    self.assertAlmostEqual(r, 0, places=16)
        
    def test1(self):
        self.do(indices=None)
    
    def test2(self):
        self.do(indices=[1])

class TestRotate(unittest.TestCase):
    def do(self, indices=None, stepsize=.8):
        natoms = 4
        x = np.array([rotations.random_aa() for _ in range(natoms)])
        x0 = x.copy()
        x = x.reshape(-1)
        
        rotate(stepsize, x, indices=indices)
        
        x = x.reshape(-1,3)
        rlist = [relative_angle(x[i], x0[i]) for i in range(natoms)]
        for r in rlist:
            self.assertLessEqual(r, stepsize)
        
        if indices is not None:
            for i, r in enumerate(rlist):
                if i in indices:
                    self.assertGreater(r, 0)
                else:
                    self.assertAlmostEqual(r, 0, places=16)
        
    def test1(self):
        self.do(indices=None)
    
    def test2(self):
        self.do(indices=[1])

    def test3(self):
        # test for small angles
        self.do(stepsize=.001)

class TestReducedCoordsDisplace(unittest.TestCase):
    def do(self, indices=None):
        # simple test with lattice matrix the identity
        natoms = 4
        x = np.zeros(natoms*3)
        
        stepsize=1.5
        scale = 0.1 # just to make things a bit more interesting
        lattice_matrix = np.eye(3) * scale
        reduced_coordinates_displace(stepsize, lattice_matrix, x, indices=indices)
        
        x = x.reshape(-1,3)
        rlist = np.array([np.linalg.norm(r) for r in x])
        rlist *= scale # undo the scaling effects
        for r in rlist:
            self.assertLessEqual(r, stepsize)
        
        if indices is not None:
            for i, r in enumerate(rlist):
                if i in indices:
                    self.assertGreater(r, 0)
                else:
                    self.assertAlmostEqual(r, 0, places=16)
        
    def test1(self):
        self.do(indices=None)
    
    def test2(self):
        self.do(indices=[1])

        
        

if __name__ == "__main__":
    unittest.main()

