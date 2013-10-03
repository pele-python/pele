import unittest
import numpy as np

from pele.systems import LJCluster
from pele.transition_states._find_lowest_eig import FindLowestEigenVector, analyticalLowestEigenvalue, findLowestEigenVector

class TestFindLowestEigenvector(unittest.TestCase):
    def setUp(self):
        self.setUp1()

    def setUp1(self, **kwargs):
        natoms = 18
#        s = LJCluster(natoms)
#        nfrozen = 6
#        reference_coords = s.get_random_configuration()
#        self.system = LJClusterFrozen(13, range(nfrozen), reference_coords)
        self.system = LJCluster(natoms)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()
        
        self.finder = FindLowestEigenVector(self.x.copy(), self.pot, **kwargs)
        

    def test(self):
        lval, lvec = analyticalLowestEigenvalue(self.x, self.pot)
        ret = self.finder.run(100)
        self.assertLess(np.abs(ret.eigenval - lval) / np.abs(lval), 1e-2)
    
        
    def test2(self):
        lval, lvec = analyticalLowestEigenvalue(self.x, self.pot)
        ret = findLowestEigenVector(self.x.copy(), self.pot)
        self.assertLess(np.abs(ret.eigenval - lval) / np.abs(lval), 1e-2)

class TestHEF_InvertedGradient(TestFindLowestEigenvector):
    def setUp(self):
        self.setUp1(inverted_gradient=True)
  
    
        
if __name__ == "__main__":
    unittest.main()