import unittest
import os

import numpy as np

from pele.thermodynamics._normalmodes import logproduct_freq2, normalmodes
from pele.thermodynamics import get_thermodynamic_information
from pele.systems import LJCluster
from pele.thermodynamics.heat_capacity import minima_to_cv, dos_to_cv


class TestHeatCapacity(unittest.TestCase):
    def setUp(self):
        import pele.rates.tests.__init__ as f

        dirname = os.path.dirname(f.__file__)
        dbfname = os.path.join(dirname, "lj15.sqlite")
        if not os.path.exists(dbfname):
            raise IOError("database file %s does not exist" % dbfname)
        self.system = LJCluster(15)
        self.db = self.system.create_database(dbfname, createdb=False)
        get_thermodynamic_information(self.system, self.db, nproc=None)

    def test_cv(self):
        kT = np.array([.1, 1.])
        k = self.system.get_ndof()
        minima = self.db.minima()[:10]
        cvdata = minima_to_cv(minima, kT, k)
        self.assertAlmostEqual(cvdata.U[0], -50.37139067, places=5)
        self.assertAlmostEqual(cvdata.U[1], -31.26584667, places=5)
        self.assertAlmostEqual(cvdata.Cv[0], 39.1151649, places=5)
        self.assertAlmostEqual(cvdata.Cv[1], 39.20781194, places=5)


_ldos = np.array([[-130.4352, 0.],
                  [-130.3616, 1.30645187],
                  [-130.288, 2.33024655],
                  [-130.2144, 3.48058723],
                  [-130.1408, 4.58681004],
                  [-130.0672, 5.64071519],
                  [-129.9936, 6.71946827],
                  [-129.92, 7.68181502],
                  [-129.8464, 8.67374929],
                  [-129.7728, 9.63850382],
                  [-129.6992, 10.57395044],
                  [-129.6256, 11.48807153],
                  [-129.552, 12.39595436],
                  [-129.4784, 13.28501673],
                  [-129.4048, 14.15311187],
                  [-129.3312, 15.00756998],
                  [-129.2576, 15.85543736],
                  [-129.184, 16.66813355],
                  [-129.1104, 17.47719029],
                  [-129.0368, 18.26983701]])


class TestDosToCv(unittest.TestCase):
    def test(self):
        ldos = _ldos
        k = 1
        T = np.array([.1, 1.])
        cvdata = dos_to_cv(ldos[:, 0], ldos[:, 1], T, K=k)
        self.assertAlmostEqual(cvdata.U[0], -129.39037151, places=5)
        self.assertAlmostEqual(cvdata.U[1], -129.10328581, places=5)
        self.assertAlmostEqual(cvdata.Cv[0], 9.26942529, places=5)
        self.assertAlmostEqual(cvdata.Cv[1], 0.50886517, places=5)


if __name__ == "__main__":
    unittest.main()

