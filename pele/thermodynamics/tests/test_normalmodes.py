import unittest
import os

import numpy as np

from pele.thermodynamics._normalmodes import logproduct_freq2, normalmodes
from pele.thermodynamics import get_thermodynamic_information
from pele.systems import LJCluster

class TestNormalModes(unittest.TestCase):
    def setUp(self):
        import pele.rates.tests.__init__ as f
        dirname = os.path.dirname(f.__file__)
        dbfname = os.path.join(dirname, "lj15.sqlite")
        if not os.path.exists(dbfname):
            raise IOError("database file %s does not exist" % dbfname)
        self.system = LJCluster(15)
        self.db = self.system.create_database(dbfname, createdb=False)
    
    def check(self, fvib_expected, coords, nzero, nnegative, metric=None): 
        pot = self.system.get_potential()
        hess = pot.getHessian(coords)
        freqs, modes = normalmodes(hess, metric=metric)
#         print v
        n, fvib = logproduct_freq2(freqs, nzero=nzero, nnegative=nnegative)
        self.assertAlmostEqual(fvib, fvib_expected, 4)
        self.assertEqual(n, len(coords) - nzero - nnegative)

    def test_minimum(self):
        m = self.db.minima()[0]
        self.check(m.fvib, m.coords, 6, 0)
    
    def test_transition_state(self):
        ts = self.db.transition_states()[0]
        self.check(ts.fvib, ts.coords, 6, 1)
        
    def test_metric_tensor(self):
        m = self.db.minima()[0]
        mt = np.eye(m.coords.size)
        self.check(m.fvib, m.coords, 6, 0, metric=mt)


    def test_get_thermo_info(self):
        newdb = self.system.create_database()
        new2old = dict()
        for ts in self.db.transition_states()[:5]:
            m1 = newdb.addMinimum(ts.minimum1.energy, ts.minimum1.coords)
            m2 = newdb.addMinimum(ts.minimum2.energy, ts.minimum2.coords)
            newts = newdb.addTransitionState(ts.energy, ts.coords, m1, m2)
            new2old[m1] = ts.minimum1
            new2old[m2] = ts.minimum2
            new2old[newts] = ts
        
        get_thermodynamic_information(self.system, newdb, nproc=2)
        
        for new in newdb.minima() + newdb.transition_states():
            old = new2old[new]
            self.assertAlmostEqual(new.energy, old.energy, 4)
            self.assertAlmostEqual(new.fvib, old.fvib, 4)
            self.assertEqual(new.pgorder, old.pgorder)


if __name__ == "__main__":
    unittest.main()
