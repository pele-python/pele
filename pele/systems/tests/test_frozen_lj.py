import unittest

from pele.systems.ljcluster_frozen import LJClusterFrozen
from pele.systems import LJCluster
from pele.thermodynamics import get_thermodynamic_information


class TestLJCluster(unittest.TestCase):
    def setUp(self):
        self.prepare_system(frozen_atoms=[0, 2, 4, 6])

    def prepare_system(self, frozen_atoms=None):
        if not frozen_atoms: frozen_atoms = [0, 2, 4]
        self.natoms = 13

        fsys = LJCluster(self.natoms)

        self.reference_coords = fsys.get_random_configuration()

        self.system = LJClusterFrozen(self.natoms, frozen_atoms, self.reference_coords)

    def make_database(self, nminima=2):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        while db.number_of_minima() < 2:
            bh.run(1)

        return db

    def test1(self):
        coords = self.system.coords_converter.get_reduced_coords(self.reference_coords)
        pot = self.system.get_potential()
        pot.getEnergy(coords)

    def test_mobile(self):
        self.assertEqual(3 * self.system.nmobile, self.system.coords_converter.get_mobile_dof().size)

    def testpermlist(self):
        permlist = self.system.get_permlist()
        self.assertEqual(len(permlist[0]), self.natoms - len(self.system.frozen_atoms))

    def test_basinhopping(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        while db.number_of_minima() < 2:
            bh.run(1)

    def test_connect(self):
        db = self.make_database()
        min1, min2 = db.minima()[:2]
        connect = self.system.get_double_ended_connect(min1, min2, db)
        connect.connect()

    # self.assertTrue(connect.success())

    def test_fvib(self):
        db = self.make_database()
        get_thermodynamic_information(self.system, db, nproc=1)

# class TestLJCluster2(TestLJCluster):
# def setUp(self):
# self.prepare_system(frozen_atoms=[0,2])
# class TestLJCluster1(TestLJCluster):
# def setUp(self):
# self.prepare_system(frozen_atoms=[0])


if __name__ == "__main__":
    unittest.main() 
