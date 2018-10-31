import unittest

from pele.systems import BLJBulk


class TestLJClusterSystem(unittest.TestCase):
    def setUp(self):
        self.natoms = 13
        boxvec = [7.] * 3
        self.system = BLJBulk(self.natoms, boxvec)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), self.natoms)
        
        self.assertIsNotNone(db.get_property("ntypeA"))
        self.assertIsNotNone(db.get_property("boxvec"))

    def test_permlist(self):
        permlist = self.system.get_permlist()
        self.assertEqual(len(permlist), 2)

    def test_bh(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        bh.run(3)
        self.assertGreater(db.number_of_minima(), 0)


if __name__ == "__main__":
    unittest.main()
