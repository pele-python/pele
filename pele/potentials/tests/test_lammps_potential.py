import unittest

from pele.optimize import lbfgs_cpp
from pele.potentials import have_lammps


@unittest.skipUnless(have_lammps, "requires LAMMPS")
class TestLammpsPotential(unittest.TestCase):
    """the potential pele.potentials exports: compiled if it was built"""

    def get_potential_class(self):
        from pele.potentials import LAMMPSPotential

        return LAMMPSPotential

    def setUp(self):
        from lammps import lammps
        from pele.potentials import LJ

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
        self.potential = self.get_potential_class()(lmp)
        self.lj = LJ()
        self.coords = self.lmp.numpy.extract_atom('x').flatten().copy()

    def tearDown(self):
        self.lmp.close()

    def test_energy(self):
        self.assertAlmostEqual(
            self.potential.getEnergy(self.coords),
            self.lj.getEnergy(self.coords),
            places=10,
        )

    def test_minimize(self):
        ret = lbfgs_cpp(self.coords, self.potential)
        self.assertTrue(ret.success)
        self.assertAlmostEqual(
            self.potential.getEnergy(ret.coords),
            self.lj.getEnergy(ret.coords),
            places=10,
        )


class TestLammpsPotentialPython(TestLammpsPotential):
    """the pure Python fallback, used when the compiled potential is not built"""

    def get_potential_class(self):
        from pele.potentials.lammps_potential_python import LAMMPSPotential

        return LAMMPSPotential
