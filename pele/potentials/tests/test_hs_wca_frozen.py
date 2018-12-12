from __future__ import division

import unittest
import numpy as np

from pele.potentials import HS_WCA
from pele.optimize import ModifiedFireCPP

class Test2dMinimization(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.L_mobile = 4
        self.L_total = self.L_mobile + 2
        self.nr_particles_mobile = self.L_mobile * self.L_mobile
        self.nr_particles_total = self.L_total * self.L_total
        self.nr_particles_frozen = self.nr_particles_total - self.nr_particles_mobile
        self.box_dimension = 2
        self.ndof = self.nr_particles_total * self.box_dimension
        self.n_frozen_dof = self.nr_particles_frozen * self.box_dimension
        self.frozen_dof = []
        self.frozen_atoms = []
        for particle_index in range(self.nr_particles_total):
            xmean = int(particle_index % self.L_total)
            ymean = int(particle_index / self.L_total)
            if ymean == 0 or ymean == self.L_total - 1 or xmean == 0 or xmean == self.L_total - 1:
                self.frozen_dof.append(particle_index * self.box_dimension)
                self.frozen_dof.append(particle_index * self.box_dimension + 1)
                self.frozen_atoms.append(particle_index)
        self.eps = 1
        self.x = np.zeros(self.ndof)
        for p in range(self.nr_particles_total):
            xmean = int(p % self.L_total)
            ymean = int(p / self.L_total)
            self.x[p * self.box_dimension] = xmean + 0.1 * np.random.rand()
            self.x[p * self.box_dimension + 1] = ymean + 0.1 * np.random.rand()
        self.radii = np.asarray([0.3 + 0.01 * np.random.rand() for _ in range(self.nr_particles_total)])
        self.sca = 1
        self.rcut = 2 * (1 + self.sca) * np.amax(self.radii)
        self.boxvec = (self.L_total + self.rcut) * np.ones(self.box_dimension)
        self.pot_cells_N_frozen_N = HS_WCA(eps=self.eps, sca=self.sca,
                                    radii=self.radii, ndim=self.box_dimension,
                                    boxvec=self.boxvec, use_periodic=True,
                                    use_frozen=False, use_cell_lists=False)
        self.pot_cells_Y_frozen_N = HS_WCA(eps=self.eps, sca=self.sca,
                                    radii=self.radii, ndim=self.box_dimension,
                                    boxvec=self.boxvec, use_periodic=True,
                                    use_frozen=False, use_cell_lists=True,
                                    reference_coords=self.x, rcut=self.rcut)
        self.pot_cells_N_frozen_Y = HS_WCA(eps=self.eps, sca=self.sca,
                                    radii=self.radii, ndim=self.box_dimension,
                                    boxvec=self.boxvec, use_periodic=True,
                                    use_frozen=True, use_cell_lists=False,
                                    frozen_atoms=self.frozen_atoms,
                                    reference_coords=self.x)
        self.pot_cells_Y_frozen_Y = HS_WCA(eps=self.eps, sca=self.sca,
                                    radii=self.radii, ndim=self.box_dimension,
                                    boxvec=self.boxvec, use_periodic=True,
                                    use_frozen=True, use_cell_lists=True,
                                    reference_coords=self.x,
                                    frozen_atoms=self.frozen_atoms, rcut=self.rcut)
        self.x_red = []
        for atom in range(self.nr_particles_total):
            if atom not in self.frozen_atoms:
                self.x_red.extend(self.x[atom * self.box_dimension : (atom + 1) * self.box_dimension])
        self.opt_NN = ModifiedFireCPP(self.x, self.pot_cells_N_frozen_N)
        self.opt_YN = ModifiedFireCPP(self.x, self.pot_cells_Y_frozen_N)
        self.opt_NY = ModifiedFireCPP(self.x_red, self.pot_cells_N_frozen_Y)
        self.opt_YY = ModifiedFireCPP(self.x_red, self.pot_cells_Y_frozen_Y)
    def test_energies(self):
        self.res_e_before_cells_N_frozen_N = self.opt_NN.get_result()
        self.res_e_before_cells_Y_frozen_N = self.opt_YN.get_result()
        self.res_e_before_cells_N_frozen_Y = self.opt_NY.get_result()
        self.res_e_before_cells_Y_frozen_Y = self.opt_YY.get_result()
        self.assertAlmostEqual(self.res_e_before_cells_N_frozen_N.energy, self.res_e_before_cells_N_frozen_Y.energy, places=8)
        self.assertAlmostEqual(self.res_e_before_cells_Y_frozen_N.energy, self.res_e_before_cells_N_frozen_N.energy, places=8)
        self.assertAlmostEqual(self.res_e_before_cells_N_frozen_N.energy, self.res_e_before_cells_Y_frozen_Y.energy, places=8)
        self.assertAlmostEqual(self.pot_cells_N_frozen_N.getEnergy(self.x), self.pot_cells_Y_frozen_N.getEnergy(self.x), places=8)
    def test_minimization(self):
        self.opt_NN.run()
        self.opt_YN.run()
        self.opt_NY.run()
        self.opt_YY.run()
        self.res_NN = self.opt_NN.get_result()
        self.res_YN = self.opt_YN.get_result()
        self.res_NY = self.opt_NY.get_result()
        self.res_YY = self.opt_YY.get_result()
        self.assertTrue(self.res_NN.success)
        self.assertTrue(self.res_YN.success)
        self.assertTrue(self.res_NY.success)
        self.assertTrue(self.res_YY.success)
        self.assertAlmostEqual(self.res_NY.energy, self.res_YY.energy, delta=1e-10)
        self.assertAlmostEqual(self.pot_cells_N_frozen_N.getEnergy(self.res_NN.coords), self.pot_cells_Y_frozen_N.getEnergy(self.res_YN.coords), delta=1e-10)
        
if __name__ == "__main__":
    unittest.main()

