import unittest
import numpy as np

from pele.systems import LJCluster
from pele.utils.frozen_atoms import FrozenCoordsConverter, FrozenPotWrapper

class TestFrozenCoordsConverter(unittest.TestCase):
    def setUp(self):
        self.ndof = 97
        
        self.reference_coords = np.random.rand(self.ndof)
        self.frozen_dof = np.array([0, 1, 2, 10, 40])
        
        
        self.converter = FrozenCoordsConverter(self.reference_coords, self.frozen_dof)
        
    def test1(self):
        reduced_coords = self.converter.get_reduced_coords(self.reference_coords)
        self.assertEqual(self.ndof - len(self.frozen_dof), len(reduced_coords))
        
    def test2(self):
        full_coords = self.converter.get_full_coords(self.converter.get_reduced_coords(self.reference_coords))
        self.assertTrue((full_coords == self.reference_coords).all())
    

if __name__ == "__main__":
    unittest.main() 