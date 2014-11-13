import unittest
import numpy as np
import os
import logging

from pele.potentials import _lj_cpp
from pele.utils.xyz import read_xyz
import _base_test


class TestLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        self.pot = _lj_cpp.LJ()
        self.natoms = 13
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_lj13_gmin.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin = float(xyz.title)


class TestErrorPotential(unittest.TestCase):
    def setUp(self):
        self.pot = _lj_cpp._ErrorPotential()
        self.x = np.random.uniform(-1, 1, [9])

    def test(self):
        with self.assertRaises(RuntimeError):
            self.pot.getEnergy(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.getEnergyGradient(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.NumericalDerivative(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.NumericalHessian(self.x)


# with self.assertRaises(NotImplementedError):
# pot.getEnergyGradient(_xrand)

class TestLJ_CPP_NeighborList(_base_test._BaseTest):
    def setUp(self):
        np.random.seed(0)
        self.natoms = 13
        nlist = [[i, j] for i in xrange(self.natoms) for j in xrange(i + 1, self.natoms)]
        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
        self.pot = _lj_cpp.LJNeighborList(nlist)
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_lj13_gmin.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin = float(xyz.title)


if __name__ == "__main__":
    logging.basicConfig(filename='lj_cpp.log', level=logging.DEBUG)
    unittest.main()
