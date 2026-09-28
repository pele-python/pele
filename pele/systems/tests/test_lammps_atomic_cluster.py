import unittest

from pele.potentials import have_lammps


@unittest.skipUnless(have_lammps, "requires LAMMPS")
class TestLJClusterSystem(unittest.TestCase):
    def setUp(self):
        from lammps import lammps
        from pele.systems import AtomicClusterLAMMPS

        self.lmp = lmp = lammps(cmdargs="-screen none -log none".split())
        lmp.commands_string("""
        units lj
        atom_style atomic
        atom_modify sort 0 0.0
        boundary f f f
        lattice fcc 0.8442
        region box block -100 100 -100 100 -100 100
        create_box 2 box
        create_atoms 1 random 80 12345 NULL
        create_atoms 2 random 20 12346 NULL
        mass * 1.0
        pair_style lj/cut 100.0
        pair_coeff 1 1 1.0 1.0 100.0
        pair_coeff 1 2 1.5 0.8 100.0
        pair_coeff 2 2 0.5 0.88 100.0
        pair_modify shift yes
        thermo_modify norm no
        """)
        self.system = AtomicClusterLAMMPS(self.lmp)

    def test_database_property(self):
        db = self.system.create_database()
        p = db.get_property("natoms")
        self.assertIsNotNone(p)
        self.assertEqual(p.value(), 100)

    def test_permlist(self):
        self.assertEqual(
            self.system.get_permlist(),
            [list(range(0, 80)), list(range(80, 100))],
        )

    def test_bh(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        bh.run(3)
        self.assertGreater(db.number_of_minima(), 0)


if __name__ == "__main__":
    unittest.main()
