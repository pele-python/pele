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

    def test_basinhopping_max_n_minima(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(database=db, max_n_minima=2)
        bh.run(10)
        self.assertEqual(db.number_of_minima(), 2)

    def test_basinhopping_max_n_minima_params(self):
        db = self.system.create_database()
        self.system.params.basinhopping.max_n_minima = 2
        bh = self.system.get_basinhopping(database=db)
        bh.run(10)
        self.assertEqual(db.number_of_minima(), 2)


if __name__ == "__main__":
    unittest.main()
