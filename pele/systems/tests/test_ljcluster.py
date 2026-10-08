import unittest

import numpy as np

from pele.systems import LJCluster


class TestLJClusterSystem(unittest.TestCase):
    def setUp(self):
        self.addCleanup(np.random.set_state, np.random.get_state())
        np.random.seed(0)
        self.natoms = 13
        self.system = LJCluster(self.natoms)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), 13)

    def test_basinhopping_max_n_minima(self):
        db = self.system.create_database()
        # Test the storage cap independently of Metropolis acceptance.
        bh = self.system.get_basinhopping(
            database=db, max_n_minima=2, insert_rejected=True
        )
        bh.run(10)
        self.assertEqual(db.number_of_minima(), 2)

    def test_basinhopping_max_n_minima_params(self):
        db = self.system.create_database()
        self.system.params.basinhopping.max_n_minima = 2
        bh = self.system.get_basinhopping(database=db, insert_rejected=True)
        bh.run(10)
        self.assertEqual(db.number_of_minima(), 2)


if __name__ == "__main__":
    unittest.main()
