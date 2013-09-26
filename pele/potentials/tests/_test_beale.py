import unittest
import numpy as np
from pele.potentials.test_functions._beale import BealeSystem

import _base_test


class TestBeale(_base_test._BaseTest):
    def setUp(self):
        self.system = BealeSystem()
        self.pot = self.system.get_potential()
        self.xmin = self.pot.target_coords
        self.Emin = self.pot.target_E
        self.xrandom = self.system.get_random_configuration()
        

if __name__ == "__main__":
    unittest.main()