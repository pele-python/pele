import unittest

import numpy as np

from pele.potentials.ljcut import LJCut


class LJCutTest(unittest.TestCase):
    def setUp(self):
        self.natoms = 10
        self.coords = np.random.uniform(-1, 1., 3 * self.natoms) * self.natoms ** (-1. / 3)
        self.pot = LJCut()
        self.E = self.pot.getEnergy(self.coords)
        self.Egrad, self.grad = self.pot.getEnergyGradient(self.coords)

        self.ilist = []  # np.zeros([self.natoms*(self.natoms-1)/2, 2])
        k = 0
        for i in range(self.natoms):
            for j in range(i):
                k += 1
                self.ilist.append([i, j])
        self.ilist = np.array(self.ilist)
        # print self.ilist

    def test_grad_e(self):
        self.assertAlmostEqual(self.E, self.Egrad, 7)

    def test_lists_e(self):
        e = self.pot.getEnergyList(self.coords, self.ilist)
        self.assertAlmostEqual(self.E, e, 7)

    def test_lists_eg(self):
        e, g = self.pot.getEnergyGradientList(self.coords, self.ilist)
        self.assertAlmostEqual(self.E, e, 7)
        gdiffmax = np.max(np.abs(g - self.grad)) / np.max(np.abs(self.grad))
        self.assertLess(gdiffmax, 1e-7)


if __name__ == "__main__":
    unittest.main()

