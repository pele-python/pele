from __future__ import absolute_import
import unittest
import logging

import numpy as np

from pele.potentials import _wca_cpp
from pele.optimize._quench import lbfgs_cpp
from . import _base_test


def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms


class TestWCA_CPP(_base_test._BaseTest):
    def setUp(self):
        boxv = [3, 3, 3]
        self.pot = _wca_cpp.WCA(ndim=3, boxvec=boxv)
        self.natoms = 25
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 0.001

        xyz = np.random.uniform(-1, 1, [3 * self.natoms])
        xyz = minimize(xyz, self.pot)
        xyz = minimize(xyz[0], self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])


class TestWCA2D_CPP(_base_test._BaseTest):
    def setUp(self):
        boxv = [3, 3]
        self.pot = _wca_cpp.WCA(ndim=2, boxvec=boxv)
        self.natoms = 25
        self.xrandom = np.random.uniform(-1, 1, [2 * self.natoms]) * 0.01
        xy = np.random.uniform(-1, 1, [2 * self.natoms])
        xy = minimize(xy, self.pot)
        xy = minimize(xy[0], self.pot)
        self.xmin = xy[0].reshape(-1).copy()
        self.Emin = float(xy[1])


class TestWCA_CPP_NeighborList(_base_test._BaseTest):
    def setUp(self):
        self.natoms = 13
        nlist = [[i, j] for i in range(self.natoms) for j in range(i + 1, self.natoms)]
        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
        self.pot = _wca_cpp.WCANeighborList(nlist)
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        xyz = minimize(self.xrandom, self.pot)
        self.xmin = xyz[0].reshape(-1).copy()
        self.Emin = float(xyz[1])


class TestWCA_CPP_AtomList(_base_test._TestConfiguration):
    def setUp(self):
        self.natoms = 13
        atoms = np.array(list(range(self.natoms - 1)))
        self.pot = _wca_cpp.WCAAtomList(atoms)
        self.x0 = np.array([-0.4225717, 1.99835681, -0.76517552, 0.59082277, -0.84699246,
                            -0.76256685, 1.57673137, -1.21902228, 1.27521222, -0.50325388,
                            -1.69997349, 0.65583764, 1.32327954, -1.81756898, 1.57247512,
                            -1.03150491, 1.41521019, 0.55743455, 1.66775241, 1.60385959,
                            -0.50645429, 0.71341477, -0.41636407, 1.36314406, -0.39649335,
                            -0.64088725, 0.27695302, -0.36016137, 1.29213068, 0.92494101,
                            0.37140092, -1.61146783, 1.75448354, 0.96222954, -0.06410995,
                            -0.32505948, 1.21724737, -1.56051696, -1.36116059])
        self.e0 = 266.7260923712355

    def test_grad_is_zero(self):
        e, g = self.pot.getEnergyGradient(self.x0)
        g = g.reshape([-1, 3])
        self.assertAlmostEqual(g[-1, 0], 0)
        self.assertAlmostEqual(g[-1, 1], 0)
        self.assertAlmostEqual(g[-1, 2], 0)


if __name__ == "__main__":
    logging.basicConfig(filename='lj_cpp.log', level=logging.DEBUG)
    unittest.main()

