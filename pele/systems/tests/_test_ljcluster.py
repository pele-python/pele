import unittest

from pele.systems import LJCluster

class TestLJClusterSystem(unittest.TestCase):
    def setUp(self):
        self.natoms = 13
        self.system = LJCluster(self.natoms)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), 13)

if __name__ == "__main__":
    unittest.main()