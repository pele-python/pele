import unittest

from pele.systems import BLJCluster


class TestBLJClusterSystem(unittest.TestCase):
    def setUp(self):
        self.natoms = 13
        self.ntypeA = 8
        self.system = BLJCluster(self.natoms, self.ntypeA)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), self.natoms)
        p = db.get_property("ntypeA")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), self.ntypeA)

# for name in ["sigAA", "sigBB", "sigAB", "epsAA", "epsAB", "epsBB"]:
#             p = db.get_property("ntypeA")
#             self.assertIsNotNone(p)

if __name__ == "__main__":
    unittest.main()
