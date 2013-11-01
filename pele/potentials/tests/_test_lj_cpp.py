import unittest
import numpy as np
import os

from pele.potentials._lj_cpp import LJ
from pele.utils.xyz import read_xyz
import _base_test

class TestLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        self.pot = LJ() 
        self.natoms = 13
        self.xrandom = np.random.uniform(-1,1,[3*self.natoms]) *5.
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_lj13_gmin.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin = float(xyz.title)

if __name__ == "__main__":
    unittest.main()
