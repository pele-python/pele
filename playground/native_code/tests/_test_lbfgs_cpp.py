import unittest
import numpy as np

from playground.native_code._lbfgs import LBFGS_CPP
from playground.native_code import _lj

class TestLBFGS_CPP(unittest.TestCase):
    def setUp(self):
        self.natoms = 18
        self.pot = _lj.LJ() 
        self.x0 = np.random.uniform(-1,1, 3*self.natoms)
        self.lbfgs = LBFGS_CPP(self.pot, self.x0)
        
    def test(self):
        ret = self.lbfgs.run()
        self.assertTrue(ret.success)
        e, g = self.pot.getEnergyGradient(ret.coords)
        self.assertAlmostEqual(e, ret.energy, delta=1e-6)
        self.assertLess(np.max(np.abs(g - ret.grad)), 1e-6)
        rms = np.linalg.norm(g) / np.sqrt(g.size)
        self.assertAlmostEqual(rms, ret.rms, 1e-6)
        
        self.assertGreater(ret.nfev, 0)
        self.assertGreater(ret.nsteps, 0)

if __name__ == "__main__":
    unittest.main()