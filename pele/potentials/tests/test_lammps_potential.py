import unittest

import numpy as np

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

    def test_gradient_with_atom_sorting(self):
        """the gradient stays in input order although LAMMPS sorts atoms by default"""
        from lammps import lammps

        lmp = lammps(cmdargs="-screen none -log none".split())
        self.addCleanup(lmp.close)
        lmp.commands_string("""
        units lj
        atom_style atomic
        boundary s s s
        region box block -6 6 -6 6 -6 6
        create_box 1 box
        create_atoms 1 random 60 12345 NULL overlap 0.9
        mass 1 1.0
        pair_style lj/cut 2.5
        pair_coeff 1 1 1.0 1.0 2.5
        """)
        potential = self.get_potential_class()(lmp)
        coords = lmp.numpy.extract_atom("x").flatten().copy()
        rng = np.random.default_rng(0)
        for _ in range(3):
            x = coords + rng.normal(scale=0.1, size=coords.size)
            _, grad = potential.getEnergyGradient(x)
            numerical = potential.NumericalDerivative(x, eps=1e-6)
            np.testing.assert_allclose(grad, numerical, rtol=1e-4, atol=1e-4)


class TestLammpsPotentialPython(TestLammpsPotential):
    """the pure Python fallback, used when the compiled potential is not built"""

    def get_potential_class(self):
        from pele.potentials.lammps_potential_python import LAMMPSPotential

        return LAMMPSPotential
