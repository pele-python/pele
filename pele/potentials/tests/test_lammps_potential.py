import unittest

from pele.optimize import lbfgs_cpp
from pele.potentials import have_lammps


@unittest.skipUnless(have_lammps, "requires LAMMPS")
class TestLammpsPotential(unittest.TestCase):
    def setUp(self):
        from lammps import lammps
        from pele.potentials import LAMMPSPotential, LJ

        self.lmp = lmp = lammps(cmdargs="-screen none -log none".split())
        lmp.cmd.units("lj")
        lmp.cmd.atom_style("atomic")
        lmp.cmd.boundary('s', 's', 's')
        lmp.cmd.lattice("fcc", 0.8442)
        lmp.cmd.region("box block", 0, 1, 0, 1, 0, 1)
        lmp.cmd.create_box(1, "box")
        lmp.cmd.create_atoms(1, "random", 13, 12345, "NULL")
        lmp.cmd.mass(1, 1.0)
        lmp.cmd.pair_style("lj/cut", 100.0)
        lmp.cmd.pair_coeff(1, 1, 1.0, 1.0, 100.0)
        lmp.cmd.neighbor(0.3, "bin")
        self.potential = LAMMPSPotential(lmp)
        self.lj = LJ()
        self.coords = self.lmp.numpy.extract_atom('x').flatten().copy()

    def tearDown(self):
        self.lmp.close()

    def test_energy(self):
        self.assertEqual(
            self.potential.getEnergy(self.coords),
            self.lj.getEnergy(self.coords)
        )

    def test_minimize(self):
        ret = lbfgs_cpp(self.coords, self.potential)
        self.assertTrue(ret.success)
        self.assertEqual(
            self.potential.getEnergy(ret.coords),
            self.lj.getEnergy(ret.coords)
        )
