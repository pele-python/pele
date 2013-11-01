import unittest
import numpy as np
import os

from pele.potentials._lj_cpp import BLJCut
from pele.utils.xyz import read_xyz
import _base_test

class TestLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_blj13_min.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        ntypeA, self.Emin, rcut, epsAA, sigAA, epsBB, sigBB, epsAB, sigAB = map(float, xyz.title.split()[1::2])
        ntypeA = int(ntypeA)
        
        natoms = self.xmin.size / 3
        
        self.pot = BLJCut(natoms, ntypeA, rcut=rcut, sigAA=sigAA, epsAA=epsAA, 
                          epsBB=epsBB, sigBB=sigBB, epsAB=epsAB, sigAB=sigAB) 
        self.xrandom = np.random.uniform(-1,1,self.xmin.size) *5.
    

if __name__ == "__main__":
    unittest.main()
