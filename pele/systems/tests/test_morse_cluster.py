import unittest

from pele.systems import MorseCluster


class TestLJClusterSystem(unittest.TestCase):
    def setUp(self):
        self.natoms = 13
        self.system = MorseCluster(self.natoms)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), 13)

    def test_permlist(self):
        permlist = self.system.get_permlist()
        self.assertEqual(len(permlist), 1)
        for i, v in enumerate(sorted(permlist[0])):
            self.assertEqual(i, v)

    def test_bh(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        bh.run(3)
        self.assertGreater(db.number_of_minima(), 0)


if __name__ == "__main__":
    unittest.main()
