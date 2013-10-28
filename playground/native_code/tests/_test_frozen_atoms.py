import unittest
import numpy as np

from playground.native_code import _lj
from pele.systems import LJClusterFrozen
from pele.optimize import mylbfgs
import _test_lj_cpp

class TestFrozenAtoms(_test_lj_cpp.TestLJ_CPP):
    def setUp(self):
        self.natoms = 18
        frozen_atoms = np.array([1, 2, 5, 9], dtype=int)
        reference_coords = np.random.uniform(-1,1, 3*self.natoms)
        self.pot = _lj.LJFrozen(reference_coords, frozen_atoms)
        
        system = LJClusterFrozen(self.natoms, frozen_atoms, reference_coords)
        self.pot_comp = system.get_potential()
        x = system.get_random_configuration()
        ret = mylbfgs(x, self.pot_comp, tol=10.)
        self.x = ret.coords
        
    
    def test(self):
        eonly = self.pot.getEnergy(self.x)
        e, g = self.pot.getEnergyGradient(self.x)
        self.assertAlmostEqual(e, eonly, delta=1e-6)
        et, gt = self.pot_comp.getEnergyGradient(self.x)
        self.assertAlmostEqual(e, et, delta=1e-6)
        self.assertLess(np.max(np.abs(g - gt)), 1e-6)

if __name__ == "__main__":
    unittest.main()
