import unittest
import numpy as np

from pele.systems import LJCluster, LJClusterFrozen, ljcluster_frozen
from pele.transition_states.find_lowest_eig import LowestEigPot

class TestEigPot(unittest.TestCase):
    def setUp(self):
        natoms = 18
#        s = LJCluster(natoms)
#        nfrozen = 6
#        reference_coords = s.get_random_configuration()
#        self.system = LJClusterFrozen(13, range(nfrozen), reference_coords)
        self.system = LJCluster(natoms)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()
        self.eigpot = LowestEigPot(self.x, self.pot, orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors())
    
    def test_gradient(self):
        e1 = self.eigpot.getEnergy(self.x)
        e, g = self.eigpot.getEnergyGradient(self.x)
        gnum = self.eigpot.NumericalDerivative(self.x)
        
        print np.max(np.abs(g)), np.max(np.abs(gnum))
        self.assertLess(np.max(np.abs(g-gnum)) / np.max(np.abs(g)), 1e-6)
        self.assertAlmostEqual(e, e1, 6)
        
if __name__ == "__main__":
    unittest.main()