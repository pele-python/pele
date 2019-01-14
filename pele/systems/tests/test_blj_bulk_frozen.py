import unittest
import numpy as np

from pele.systems import BLJBulkFrozen


class TestLJClusterSystem(unittest.TestCase):
    def setUp(self):
        self.natoms = 40
        boxl = 3.
        boxvec = [boxl] * 3
        reference_coords = np.random.uniform(0,-boxl, self.natoms*3)
        frozen_atoms = [0, 1, 5, 7, self.natoms-1]
        self.system = BLJBulkFrozen(self.natoms, boxvec, reference_coords, frozen_atoms)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), self.natoms)
        
        self.assertIsNotNone(db.get_property("ntypeA"))
        self.assertIsNotNone(db.get_property("boxvec"))
        self.assertIsNotNone(db.get_property("reference_coords"))
        self.assertIsNotNone(db.get_property("frozen_atoms"))


    def test_permlist(self):
        permlist = self.system.get_permlist()
        self.assertEqual(len(permlist), 2)

    def test_bh(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        bh.run(5)
        self.assertGreater(db.number_of_minima(), 0)


if __name__ == "__main__":
    unittest.main()
