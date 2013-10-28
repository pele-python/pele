import unittest
import numpy as np

from playground.native_code import _lj
from pele.systems import LJCluster
from pele.potentials import LJpshift
from pele.optimize import mylbfgs
import _test_lj_cpp

class TestBLJ_CPP(_test_lj_cpp.TestLJ_CPP):
    def setUp(self):
        self.natoms = 18
        ntypeA = 5
        self.pot = _lj.BLJCut(self.natoms, ntypeA) 
        
        self.pot_comp = LJpshift(self.natoms, ntypeA)
        x = np.random.uniform(-1,1, 3*self.natoms)
        ret = mylbfgs(x, self.pot_comp, tol=10.)
        self.x = ret.coords
        

if __name__ == "__main__":
    unittest.main()
