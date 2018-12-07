import unittest
import os

import numpy as np

from pele.systems import LJCluster
from pele.transition_states._find_lowest_eig import LowestEigPot


class TestEigPot(unittest.TestCase):
    def setUp(self):
        self.first_order = False
        self.setUp1()
    def setUp1(self):
        np.random.seed(0)
        natoms = 18
#        s = LJCluster(natoms)
#        nfrozen = 6
#        reference_coords = s.get_random_configuration()
#        self.system = LJClusterFrozen(13, range(nfrozen), reference_coords)
        self.system = LJCluster(natoms)
        self.x = self.system.get_random_configuration()
        #ret = self.system.get_random_minimized_configuration(tol=100.)
        #self.x = ret.coords
        self.pot = self.system.get_potential()
        self.eigpot = LowestEigPot(self.x, self.pot, 
                                   orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
                                   first_order=self.first_order,
                                   )
    
    def test_gradient(self):
        vec = np.random.rand(self.x.size)
        vec /= np.linalg.norm(vec)
        e1 = self.eigpot.getEnergy(vec)
        e, g = self.eigpot.getEnergyGradient(vec)
        gnum = self.eigpot.NumericalDerivative(vec, eps=1e-6)
        gnum -= np.dot(gnum, vec)
        
        self.assertLess(np.max(np.abs(g-gnum)) / np.max(np.abs(g)), 1e-2)
        self.assertAlmostEqual(e, e1, delta=e * 1e-4)
        self.assertAlmostEqual(1., np.dot(g, gnum) / np.linalg.norm(g) / np.linalg.norm(gnum), 3)
        
    def test_ts(self):
        from pele.utils.xyz import read_xyz
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = read_xyz(open(path + "/lj18_ts.xyz", "r"))
        x = xyz.coords.flatten()
        
        vec = np.random.rand(x.size)
        vec /= np.linalg.norm(vec)
        self.eigpot = LowestEigPot(x, self.pot, 
                                   orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
                                   first_order=self.first_order)

        e1 = self.eigpot.getEnergy(vec)
        e, g = self.eigpot.getEnergyGradient(vec)
        gnum = self.eigpot.NumericalDerivative(vec, eps=1e-4)
        gnum -= np.dot(gnum, vec)
        
        self.assertLess(np.max(np.abs(g-gnum)) / np.max(np.abs(g)), 1e-3)
        self.assertAlmostEqual(e, e1, delta=e * 1e-4)
        self.assertAlmostEqual(1., np.dot(g, gnum) / np.linalg.norm(g) / np.linalg.norm(gnum), 3)

class TestEigPotFirstOrder(TestEigPot):
    def setUp(self):
        self.first_order = True
        self.setUp1()

        
if __name__ == "__main__":
    unittest.main()

