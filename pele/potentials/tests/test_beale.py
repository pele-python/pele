from __future__ import absolute_import
import unittest

from pele.potentials.test_functions import BealeSystem
from . import _base_test


class TestBeale(_base_test._BaseTest):
    def setUp(self):
        self.system = BealeSystem()
        self.pot = self.system.get_potential()
        self.xmin = self.pot.target_coords
        self.Emin = self.pot.target_E
        self.xrandom = self.system.get_random_configuration()
        
#        print self.xmin, self.pot.getEnergy(self.xmin), self.pot.getEnergy(np.array([3., 0.1]))
        

if __name__ == "__main__":
    unittest.main()

