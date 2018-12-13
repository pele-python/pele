import unittest
import os
import sys

import numpy as np

from pele.thermodynamics._normalmodes import logproduct_freq2, normalmodes,\
    NormalModeError
from pele.thermodynamics import get_thermodynamic_information
from pele.systems import LJCluster


class TestNormalModes(unittest.TestCase):
    def setUp(self):
        import numpy as np
        s = np.random.randint(1000000)
        s = 322846
        self.seed = s
        sys.stderr.write("setUp: seed {}\n".format(self.seed))
        np.random.seed(s)
        import pele.rates.tests.__init__ as f

        dirname = os.path.dirname(f.__file__)
        dbfname = os.path.join(dirname, "lj15.{}.sqlite".format(sys.version_info.major))
        if not os.path.exists(dbfname):
            raise IOError("database file %s does not exist" % dbfname)
        self.system = LJCluster(15)
        self.db = self.system.create_database(dbfname, createdb=False)

    def check(self, fvib_expected, coords, nzero, nnegative, metric=None):
        pot = self.system.get_potential()
        hess = pot.getHessian(coords)
        freqs, modes = normalmodes(hess, metric=metric)
        # print v
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
        # note, there is an intermittant error in this test
        # it causes the system to lock, and has to do with multiprocessing
        sys.stderr.write("test_get_thermo_info: seed {}\n".format(self.seed))
        newdb = self.system.create_database()
        new2old = dict()
        for ts in self.db.transition_states()[:5]:
            m1 = newdb.addMinimum(ts.minimum1.energy, ts.minimum1.coords)
            m2 = newdb.addMinimum(ts.minimum2.energy, ts.minimum2.coords)
            newts = newdb.addTransitionState(ts.energy, ts.coords, m1, m2)
            new2old[m1] = ts.minimum1
            new2old[m2] = ts.minimum2
            new2old[newts] = ts

        get_thermodynamic_information(self.system, newdb, nproc=2, verbose=True)

        for new in newdb.minima() + newdb.transition_states():
            old = new2old[new]
            self.assertAlmostEqual(new.energy, old.energy, 4)
            self.assertAlmostEqual(new.fvib, old.fvib, 4)
            self.assertEqual(new.pgorder, old.pgorder)

    def test_too_few_zero_modes(self):
        self.system.get_nzero_modes = lambda: 10
        newdb = self.system.create_database()
        for ts in self.db.transition_states()[:4]:
            m1 = newdb.addMinimum(ts.minimum1.energy, ts.minimum1.coords)
            m2 = newdb.addMinimum(ts.minimum2.energy, ts.minimum2.coords)
            newdb.addTransitionState(ts.energy, ts.coords, m1, m2)

        get_thermodynamic_information(self.system, newdb, nproc=2, verbose=True)
        for ts in newdb.transition_states():
            self.assertTrue(ts.invalid)
        for m in newdb.minima():
            self.assertTrue(m.invalid)

    def test_too_few_negative_modes(self):
        newdb = self.system.create_database()
        tslist = []
        for ts in self.db.transition_states()[:2]:
            m1 = newdb.addMinimum(ts.minimum1.energy, ts.minimum1.coords)
            m2 = newdb.addMinimum(ts.minimum2.energy, ts.minimum2.coords)
            # add a minima as a transition state
            newts = newdb.addTransitionState(m1.energy, m1.coords, m1, m2)
            tslist.append(newts)

        # with self.assertRaises(ValueError):
        get_thermodynamic_information(self.system, newdb, nproc=2, verbose=False)
        for ts in tslist:
            self.assertTrue(ts.invalid)

    def test_too_many_negative_modes(self):
        newdb = self.system.create_database()
        mlist = []
        for ts in self.db.transition_states()[:2]:
            # add a transition state as a minimum
            m1 = newdb.addMinimum(ts.energy, ts.coords)
            m2 = newdb.addMinimum(ts.minimum2.energy, ts.minimum2.coords)
            newts = newdb.addTransitionState(m1.energy, m1.coords, m1, m2)
            mlist.append(m1)

        # with self.assertRaises(ValueError):
        get_thermodynamic_information(self.system, newdb, nproc=2, verbose=False)
        for m in mlist:
            self.assertTrue(m.invalid)


if __name__ == "__main__":
    unittest.main()
# t = TestNormalModes(methodName="test_too_few_negative_modes")
# t.setUp()
#    t.test_too_few_negative_modes()

